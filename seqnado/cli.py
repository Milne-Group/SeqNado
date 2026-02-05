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
from typing import TYPE_CHECKING, Any, List, Optional

import click
import typer
from loguru import logger

from seqnado.outputs.multiomics import find_assay_config_paths

# Optional: prettier tracebacks/console with rich if available
try:
    from rich.traceback import install as _rich_tb_install

    _rich_tb_install(show_locals=False)
except Exception:
    pass

# Rich console capture may or may not be available — keep an opt-in safe reference
try:
    from rich.console import Console  # type: ignore
    from rich.text import Text  # type: ignore

    _RICH_CONSOLE = Console(force_terminal=True)
except Exception:
    _RICH_CONSOLE = None
    Text = None  # type: ignore

if TYPE_CHECKING:
    # Available to type-checkers only; does not execute at runtime.
    import pandas as pd

from seqnado import Assay
from seqnado._version import __version__


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


# ------------------------------- Utilities ---------------------------------- #


def _pkg_traversable(pkg: str):
    """
    Return the importlib.resources Traversable for a package.
    Use resources.as_file(...) when you need a filesystem Path for a particular file.
    """
    return resources.files(pkg)


def _read_json(path: Path) -> dict:
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
    return shutil.which("snakemake") is not None


def _configure_logging(verbose: bool) -> None:
    logger.remove()
    logger.add(sys.stderr, level="DEBUG" if verbose else "INFO", colorize=True)


def _find_fastqs(hints: List[str]) -> List[Path]:
    """
    Search the provided hint directories for *.fastq.gz files.
    Skip hints that don't exist and return a sorted list.
    Searches both directly in the hint directory and in subdirectories.
    """
    out: List[Path] = []
    for loc in hints:
        p = Path(loc)
        if not p.exists():
            continue
        # Search directly in the directory
        out.extend(sorted(p.glob("*.fastq.gz")))
        # Also search one level deeper (for fastqs/{assay}/ structure)
        out.extend(sorted(p.glob("*/*.fastq.gz")))
    return out


def _style_name_with_rich(name: str, style: str = "bold cyan") -> str:
    """Return `name` styled with Rich as an ANSI string, or plain name on fallback."""
    if _RICH_CONSOLE is None or Text is None:
        return name
    with _RICH_CONSOLE.capture() as cap:
        _RICH_CONSOLE.print(Text(name, style=style), end="")
    return cap.get()


# ---- Pandera schema helpers ------------------------------------------------- #
# (these functions are defensive and robust against Pandera API differences)


def _extract_candidate_defaults_from_schema(
    model: type, assay: Any
) -> dict[str, dict[str, Any]]:
    """
    Build a metadata dict for columns that we might add with defaults.

    - Ignores: r2, r1_control, r2_control, ip, control
    - Ignores: deseq2 for non-RNA assays
    - Considers a column a candidate if it has a schema default OR is nullable
    """
    try:
        schema = model.to_schema()  # Pandera DataFrameModel -> Schema
    except Exception as e:
        logger.debug("Could not convert model to schema: %s", e)
        return {}

    out: dict[str, dict[str, Any]] = {}
    to_ignore = {"r2", "r1_control", "r2_control", "ip", "control"}

    # Also ignore deseq2 and group for non-RNA assays
    from seqnado.inputs import Assay as AssayEnum
    if assay != AssayEnum.RNA:
        to_ignore.add("deseq2")
        to_ignore.add("group")

    for name, col in schema.columns.items():
        if name in to_ignore:
            continue

        default = getattr(col, "default", None)
        nullable = bool(getattr(col, "nullable", False))
        if not (default is not None or nullable):
            continue

        description = getattr(col, "description", None)
        dtype = getattr(col, "dtype", None)

        categories = None
        try:
            import pandas as pd

            if isinstance(dtype, pd.CategoricalDtype):
                categories = list(dtype.categories)
            elif getattr(dtype, "categories", None) is not None:
                categories = list(dtype.categories)
            elif str(dtype) == "category" and name == "assay":
                try:
                    categories = [a.value for a in type(assay)]
                except Exception:
                    categories = None
        except Exception:
            categories = None

        out[name] = {
            "default": default,
            "nullable": nullable,
            "description": description,
            "dtype": dtype,
            "categories": categories,
        }

    return out


def _format_col_hint(name: str, meta: dict[str, Any]) -> str:
    """
    Build a hint like: <colored name>: <description> · choices=[...]
    """
    name_colored = _style_name_with_rich(name)

    parts: List[str] = []
    if meta.get("description"):
        parts.append(str(meta["description"]))

    cats = meta.get("categories")
    if cats:
        parts.append("choices=[" + ", ".join(map(str, cats)) + "]")

    return f"{name_colored}: " + " · ".join(parts) if parts else name_colored


def _coerce_value_to_dtype(
    value: str, dtype: Any, categories: Optional[List[Any]]
) -> Any:
    """
    Best-effort conversion from string input to the schema dtype.
    Keeps string if unsure. Enforces categorical choices if provided.
    """
    import pandas as pd
    from pandas.api import types as _pdt

    if categories is not None:
        # Allow empty/nullable
        if value == "" or value in categories:
            return value
        raise ValueError(f"Value must be one of: {', '.join(map(str, categories))}")

    if dtype is None:
        return value

    try:
        pd_dtype = pd.api.types.pandas_dtype(dtype)
    except Exception:
        pd_dtype = dtype

    # boolean handling
    if _pdt.is_bool_dtype(pd_dtype):
        low = value.strip().lower()
        if low in {"true", "t", "yes", "y", "1"}:
            return True
        if low in {"false", "f", "no", "n", "0"}:
            return False
        raise ValueError("Enter a boolean (y/n, true/false, 1/0).")

    # integer
    if _pdt.is_integer_dtype(pd_dtype):
        try:
            return int(value)
        except Exception as e:
            raise ValueError(str(e))

    # float
    if _pdt.is_float_dtype(pd_dtype):
        try:
            return float(value)
        except Exception as e:
            raise ValueError(str(e))

    # fallback to string
    return value


def _extract_deseq2_groups_from_sample_names(
    sample_ids: pd.Series,
    pattern: str | None = None,
) -> pd.Series | None:
    """
    Extract DESeq2 group labels from sample names.

    Supports multiple extraction strategies:
    1. Regex pattern matching (if pattern provided) - extracts first capture group
    2. Keyword detection (control, treated, wt, ko, etc.)
    3. Heuristic-based extraction (part before "rep" pattern)

    Args:
        sample_ids: Series of sample names
        pattern: Optional regex pattern to extract groups. The first capture group will be used.
                 Example: r'-(\\w+)-rep' for "sample-GROUPNAME-rep1"
                 Example: r'^([^-]+)' for first part before hyphen

    Returns:
        Series of group labels, or None if groups cannot be reliably extracted
    """
    import re

    # Strategy 1: Regex pattern matching
    if pattern:
        try:
            groups = sample_ids.str.extract(pattern, expand=False)
            if groups.notna().all() and groups.nunique() >= 2:
                # Validate groups are reasonably balanced
                group_counts = groups.value_counts()
                if group_counts.min() > 0 and group_counts.max() / len(groups) <= 0.9:
                    return groups
        except Exception:
            pass  # Fall through to other strategies

    # Strategy 2 & 3: Keyword detection and heuristics
    # Common group keywords to look for
    group_keywords = [
        "control", "ctrl",
        "treated", "treatment", "treat",
        "wt", "wildtype", "wild-type",
        "ko", "knockout", "knock-out",
        "mut", "mutant",
        "untreated",
        "knockdown", "kd",
        "vehicle",
        "mock",
        "dmso",
    ]

    def extract_group(sample_id: str) -> str:
        """Extract group from a single sample_id."""
        sample_lower = sample_id.lower()

        # Split on both hyphens and underscores
        parts = re.split(r'[-_]', sample_id)

        # First, try to find exact keyword matches in the parts
        for keyword in group_keywords:
            for part in parts:
                if keyword == part.lower():
                    return part

        # If no exact match, look for keywords as substrings
        for keyword in group_keywords:
            if keyword in sample_lower:
                # Try to extract the part containing the keyword
                for part in parts:
                    if keyword in part.lower():
                        return part

        # Try to find the part(s) before "rep" pattern (be more specific: "rep" followed by digits)
        for i, part in enumerate(parts):
            if re.match(r"^rep\d+$", part.lower()) and i > 0:
                # Strategy: Collect all meaningful parts between start and rep
                # Skip batch numbers and generic prefixes like "test", "sample", "exp"
                group_parts = []
                found_keyword = False

                for j in range(i - 1, -1, -1):
                    prev_part = parts[j]
                    prev_lower = prev_part.lower()

                    # Skip batch numbers - they're technical, not biological groups
                    if re.match(r"^batch\d+$", prev_lower):
                        continue

                    # Check if this is a generic prefix to stop at
                    if prev_lower in ["test", "sample", "exp", "experiment"]:
                        # We've gone too far - stop here
                        break

                    # Add this part to our group
                    group_parts.insert(0, prev_part)

                    # Track if we found a treatment keyword
                    if any(keyword in prev_lower for keyword in group_keywords):
                        found_keyword = True

                # If we collected parts, join them with hyphens
                if group_parts:
                    return "-".join(group_parts)
                # Fallback: use the part right before rep
                return parts[i - 1]

        # Fallback: use second-to-last part if there are multiple parts
        if len(parts) >= 2:
            return parts[-2]

        # Last resort: return the whole sample_id (will likely fail validation)
        return sample_id

    # Extract groups for all samples
    groups = sample_ids.apply(extract_group)

    # Check if we have at least 2 distinct groups
    unique_groups = groups.nunique()
    if unique_groups < 2:
        # Not enough groups to do differential analysis
        return None

    # Check if groups are reasonably balanced (no single group dominates too much)
    # This helps avoid extracting noise like replicate numbers
    group_counts = groups.value_counts()
    if group_counts.min() == 0 or group_counts.max() / len(groups) > 0.9:
        return None

    # Check if each "group" only has one sample (likely failed extraction)
    # This happens with samples like A1, A2, B1, B2 where each is its own "group"
    if group_counts.max() == 1:
        # Every sample is in its own group - this is not useful for DESeq2
        return None

    # Determine which group is the control (reference) and which is treatment
    # Control keywords (should be coded as 0)
    control_keywords = ["control", "ctrl", "untreated", "vehicle", "mock", "dmso", "wt", "wildtype"]

    # Find which group is the control
    unique_group_names = groups.unique()
    control_group = None

    for group_name in unique_group_names:
        if any(keyword in str(group_name).lower() for keyword in control_keywords):
            control_group = group_name
            break

    # If no control keyword found, use the first group alphabetically as reference
    if control_group is None:
        control_group = sorted(unique_group_names)[0]

    # Only create binary encoding if there are exactly 2 groups
    # For 3+ groups, DESeq2 requires more complex contrasts that should be manually specified
    if len(unique_group_names) == 2:
        # Create binary encoding: 0 for control, 1 for treatment
        deseq2_binary = groups.apply(lambda x: 0 if x == control_group else 1)
    else:
        # For 3+ groups, return None for deseq2 - user must manually configure
        logger.warning(
            f"Found {len(unique_group_names)} groups: {sorted(unique_group_names)}. "
            f"Binary DESeq2 encoding only supports 2-group comparisons. "
            f"The 'group' column will be populated, but 'deseq2' will be left empty. "
            f"Please manually specify DESeq2 levels for multi-group comparisons."
        )
        deseq2_binary = None

    # Return tuple of (group_names, binary_encoding)
    return (groups, deseq2_binary)


def _apply_interactive_defaults(
    df_in: pd.DataFrame,
    candidates: dict[str, dict[str, Any]],
    interactive: bool,
    accept_all_defaults: bool,
    deseq2_pattern: str | None = None,
    assay: Any = None,
) -> pd.DataFrame:
    """
    For any candidate column not present in df_in, optionally add it.

    - Returns a new DataFrame (does not mutate the input).
    - If `accept_all_defaults` is True, auto-add only when a schema default exists.
    - If `interactive` is False, do nothing (safe for CI/batch).
    """
    import pandas as pd

    df = df_in.copy()
    missing = [c for c in candidates.keys() if c not in df.columns]
    if not missing:
        return df

    for col in missing:
        # Skip if column was already populated during this loop
        # (e.g., both group and deseq2 populated when handling group)
        if col in df.columns:
            continue

        meta = candidates[col]
        default = meta.get("default", None)
        nullable = bool(meta.get("nullable", False))
        categories = meta.get("categories")
        dtype = meta.get("dtype")

        hint = _format_col_hint(col, meta)

        # Special handling for deseq2 and group columns: try to extract groups from sample names (RNA only)
        if accept_all_defaults and col in ("deseq2", "group") and "sample_id" in df.columns and assay is not None:
            from seqnado.inputs import Assay as AssayEnum
            if assay == AssayEnum.RNA:
                result = _extract_deseq2_groups_from_sample_names(df["sample_id"], pattern=deseq2_pattern)
                if result is not None:
                    groups, deseq2_binary = result
                    df["group"] = groups
                    if deseq2_binary is not None:
                        df["deseq2"] = deseq2_binary
                        logger.info(f"Adding 'group' and 'deseq2' columns with groups: {sorted(set(groups))}")
                    else:
                        # 3+ groups case - populate group column and add empty deseq2
                        df["deseq2"] = pd.NA  # Add empty deseq2 column to prevent prompting
                        logger.info(f"Adding 'group' column with {len(set(groups))} groups: {sorted(set(groups))}")
                    continue

        if accept_all_defaults and default is not None:
            logger.info(f"Adding '{col}' with schema default={default!r}  ({hint})")
            try:
                df[col] = pd.Series([default] * len(df), index=df.index)
            except Exception:
                df[col] = default
            continue
        elif accept_all_defaults:
            continue

        if not interactive:
            continue

        # Interactive path
        # Special handling for group/deseq2: try to extract groups first (RNA only)
        extracted_result = None
        if col in ("group", "deseq2") and "sample_id" in df.columns and assay is not None:
            from seqnado.inputs import Assay as AssayEnum
            if assay == AssayEnum.RNA:
                extracted_result = _extract_deseq2_groups_from_sample_names(df["sample_id"], pattern=deseq2_pattern)
                if extracted_result is not None:
                    groups, deseq2_binary = extracted_result
                    unique_groups = sorted(set(groups))
                    if deseq2_binary is not None:
                        add_col = typer.confirm(
                            f"⚠️ Column '{col}' is missing.\n{hint}\n"
                            f"Auto-detected groups from sample names: {unique_groups}\n"
                            f"Use these groups or enter 'n' to specify manually?",
                            default=True,
                        )
                        if add_col:
                            df["group"] = groups
                            df["deseq2"] = deseq2_binary
                            logger.info(f"Added 'group' and 'deseq2' columns with auto-detected groups: {unique_groups}")
                            continue
                    else:
                        # 3+ groups - offer uniform value or use detected groups
                        typer.echo(f"⚠️ Auto-detected {len(unique_groups)} groups: {unique_groups}")
                        typer.echo("Multi-group comparisons require manual DESeq2 configuration.")
                        use_detected = typer.confirm(
                            f"Use auto-detected groups for 'group' column (leave 'deseq2' empty for manual config)?",
                            default=True,
                        )
                        if use_detected:
                            df["group"] = groups
                            df["deseq2"] = pd.NA  # Add empty deseq2 column to prevent prompting
                            logger.info(f"Added 'group' column with {len(unique_groups)} groups")
                            continue
                        else:
                            # Fall through to let user specify uniform value
                            pass

        if default is not None:
            add_col = typer.confirm(
                f"⚠️ Column '{col}' is missing.\n{hint}\nAdd it with default or enter 'n' to skip?",
                default=True,
            )
            if add_col:
                try:
                    df[col] = pd.Series([default] * len(df), index=df.index)
                except Exception:
                    df[col] = default
        else:
            add_col = typer.confirm(
                f"⚠️ Optional column '{col}' is missing.\n{hint}\nAdd it with a uniform value or enter 'n' to skip?",
                default=False,
            )
            if add_col:
                while True:
                    raw = typer.prompt(
                        f"Enter a default value for '{col}' "
                        f"(press Enter for empty{f'; choices: {categories}' if categories else ''})",
                        default="",
                    )
                    if raw == "" and nullable:
                        df[col] = pd.Series([pd.NA] * len(df), index=df.index)
                        break
                    try:
                        coerced = _coerce_value_to_dtype(raw, dtype, categories)
                        try:
                            df[col] = pd.Series([coerced] * len(df), index=df.index)
                        except Exception:
                            df[col] = coerced
                        break
                    except ValueError as e:
                        typer.echo(f"[invalid] {e}")

    return df


# ------------------------------ Autocomplete -------------------------------- #


def _assay_names() -> List[str]:
    # Import lazily to keep CLI startup snappy
    from seqnado.inputs import Assay  # local import

    return list(Assay.all_assay_clean_names())


def assay_autocomplete(_: str) -> List[str]:
    return _assay_names()


def fastq_autocomplete(incomplete: str) -> List[str]:
    # Provide some file suggestions without expensive globbing
    p = Path(incomplete or ".")
    base = p.parent if p.name else p
    pattern = p.name or "*"
    try:
        candidates = sorted(base.glob(pattern + "*.fastq.gz"))
    except Exception:
        candidates = []
    return [str(x) for x in candidates if x.is_file()][:40]


def _get_profile_name(fn: Path) -> Optional[str]:
    name = fn.name
    if not name.startswith("profile_"):
        return None

    profile_parts = name.split("_")
    if len(profile_parts) < 2:
        return None
    initials = "".join(part[0] for part in profile_parts[1:] if part)
    return initials


def _preset_profiles() -> List[str]:
    profiles = [
        f.name
        for f in _pkg_traversable("seqnado.workflow.envs.profiles").iterdir()
        if f.is_dir()
    ]

    # At the moment these referenced by two unique letters based on the name after profile_ e.g. 'slurm_sigularity' -> 'ss'
    profiles = {
        _get_profile_name(Path(p)): p for p in profiles if p.startswith("profile_")
    }
    return profiles


def _profile_autocomplete() -> List[str]:
    return list(_preset_profiles().keys())


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
    dry_run: bool = typer.Option(
        False, help="Show actions without writing files or running scripts."
    ),
    verbose: bool = typer.Option(
        False, "--verbose", "-v", help="Increase logging verbosity."
    ),
) -> None:
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
    profiles = _preset_profiles()
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


# -------------------------------- utils ----------------------------------- #


@app.command(
    help="Manage genome configurations (list, edit, build, or generate fastq-screen config)"
)
def genomes(
    subcommand: str = typer.Argument(..., help="list, edit, build, or fastqscreen"),
    assay: str = typer.Argument(
        "atac",
        metavar="ASSAY",
        autocompletion=assay_autocomplete,
        show_choices=True,
        help="Assay type. Options: " + ", ".join(_assay_names()),
    ),
    fasta: Optional[Path] = typer.Option(
        None, "--fasta", "-f", help="Input FASTA (required for build)"
    ),
    name: Optional[str] = typer.Option(
        None, "--name", "-n", help="Genome name (prefix) for built genome"
    ),
    outdir: Path = typer.Option(
        Path.cwd() / "genome_build", "--outdir", "-o", help="Output directory for build"
    ),
    output: Optional[Path] = typer.Option(
        None,
        "-s",
        "--screen",
        help="Output path for fastqscreen config (fastqscreen subcommand)",
    ),
    threads: int = typer.Option(
        8,
        "-t",
        "--threads",
        help="Number of threads for Bowtie2 (fastqscreen subcommand)",
    ),
    no_contaminants: bool = typer.Option(
        False,
        "--no-contaminants",
        help="Exclude contaminant databases (fastqscreen subcommand)",
    ),
    contaminant_path: Optional[Path] = typer.Option(
        None,
        "--contaminant-path",
        help="Path to contaminant reference files (fastqscreen subcommand)",
    ),
    verbose: bool = typer.Option(
        False, "--verbose", "-v", help="Increase logging verbosity"
    ),
) -> None:
    """
    Manage genome configurations.

    Subcommands:
      list        Show packaged and user genome presets
      edit        Open user genome config in $EDITOR
      build       Build a genome from FASTA (not yet implemented)
      fastqscreen Generate FastqScreen configuration file
    """
    # Import locally for snappy startup
    from seqnado.config import load_genome_configs  # local import

    sub = subcommand.lower().strip()
    if sub == "list":
        try:
            cfg = load_genome_configs(assay=assay)
        except Exception as e:
            logger.error(f"Failed to load genome config: {e}")
            raise typer.Exit(code=1)

        if not cfg:
            logger.warning("No genome config found.")
            raise typer.Exit(code=0)

        for name, details in cfg.items():
            typer.echo(f"[bold]{name}[/bold]")
            try:
                items = details.dict()
            except Exception:
                items = dict(details)
            for k, v in items.items():
                typer.echo(f"  {k}: {v or '[not set]'}")
            typer.echo("")
        raise typer.Exit(code=0)

    elif sub == "edit":
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

    elif sub == "build":
        if not fasta:
            logger.error("The --fasta option is required for 'build'.")
            raise typer.Exit(code=2)

        fasta = Path(fasta).expanduser().resolve()
        if not fasta.exists():
            logger.error(f"FASTA not found: {fasta}")
            raise typer.Exit(code=3)

        outdir = Path(outdir).expanduser().resolve()
        outdir.mkdir(parents=True, exist_ok=True)

        raise NotImplementedError("Genome build not yet implemented.")

    elif sub == "fastqscreen":
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

    else:
        logger.error("Unknown subcommand for genomes: %s", sub)
        raise typer.Exit(code=2)


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
    verbose: bool = typer.Option(
        False, "--verbose", "-v", help="Increase logging verbosity."
    ),
    interactive: bool = typer.Option(
        True,
        "--interactive/--no-interactive",
        help="Interactively prompt for config values. Non-interactive mode only works for single assay configs (except MCC and multiomics).",
    ),
) -> None:
    _configure_logging(verbose)

    # Local imports to keep help fast
    from importlib.metadata import version as _pkg_version

    from seqnado.config import (
        build_default_workflow_config,
        build_workflow_config,
        render_config,
        build_multiomics_config,
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
    if assay not in Assay.all_assay_clean_names():
        allowed = ", ".join(Assay.all_assay_clean_names())
        logger.error(f"Unknown assay '{assay}'. Allowed: {allowed} or 'multiomics'")
        raise typer.Exit(code=2)

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


# -------------------------------- design ------------------------------------ #


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
        + ", ".join(Assay.all_assay_clean_names())
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
    verbose: bool = typer.Option(
        False, "--verbose", "-v", help="Increase logging verbosity."
    ),
) -> None:
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

            _assay = AssayEnum.from_clean_name(assay_name)

            # Create design object
            if _assay in {AssayEnum.CHIP, AssayEnum.CAT}:
                design_obj = FastqCollectionForIP.from_fastq_files(
                    assay=_assay,
                    files=fastq_files,
                    ip_to_control_map=ip_to_control_map,
                )
            else:
                design_obj = FastqCollection.from_fastq_files(
                    assay=_assay, files=fastq_files
                )

            df = design_obj.to_dataframe().sort_values("sample_id")

            schema_candidates = _extract_candidate_defaults_from_schema(
                DesignDataFrame, _assay
            )
            df = _apply_interactive_defaults(
                df,
                schema_candidates,
                interactive=interactive,
                accept_all_defaults=accept_all_defaults,
                deseq2_pattern=deseq2_pattern,
                assay=_assay,
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
    elif assay in Assay.non_multiomics_assays():
        allowed = ", ".join(AssayEnum.all_assay_clean_names())
        logger.error(f"Unknown assay '{assay}'. Allowed: {allowed} or 'multiomics'")
        raise typer.Exit(code=2)

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

    _assay = AssayEnum.from_clean_name(assay)
    logger.info(f"Found {len(fastq_paths)} FASTQ files for assay '{_assay.value}'")

    # Set default output filename if not provided
    if output is None:
        output = Path(f"metadata_{_assay.clean_name}.csv")

    if _assay in {AssayEnum.CHIP, AssayEnum.CAT}:
        design_obj = FastqCollectionForIP.from_fastq_files(
            assay=_assay, 
            files=fastq_paths,
            ip_to_control_map=ip_to_control_map,
        )
    else:
        design_obj = FastqCollection.from_fastq_files(assay=_assay, files=fastq_paths)

    df = design_obj.to_dataframe().sort_values("sample_id")

    schema_candidates = _extract_candidate_defaults_from_schema(DesignDataFrame, _assay)
    df = _apply_interactive_defaults(
        df,
        schema_candidates,
        interactive=interactive,
        accept_all_defaults=accept_all_defaults,
        deseq2_pattern=deseq2_pattern,
        assay=_assay,
    )

    if group_by:
        if isinstance(group_by, str) and group_by in df.columns:
            df["consensus_group"] = df[group_by].astype(str)
            logger.info(
                f"Grouped samples by column '{group_by}' into 'consensus_group'."
            )
        else:
            # Treat group_by as a regex pattern to extract from sample_id
            try:
                if _assay in [AssayEnum.CHIP, AssayEnum.CAT]:
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


# -------------------------------- pipeline ---------------------------------- #
# Allow pass-through of *unknown* options to Snakemake via ctx.args
# Replace your current pipeline function with the following:

# Allow pass-through of *unknown* options to Snakemake via ctx.args
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
    preset: str = typer.Option(
        "le",
        "--preset",
        click_type=click.Choice(_profile_autocomplete(), case_sensitive=False),
        help="Snakemake job profile preset.",
        case_sensitive=False,
    ),
    clean_symlinks: bool = typer.Option(
        False, help="Remove symlinks created by previous runs."
    ),
    scale_resources: float = typer.Option(
        1.0, "-s", "--scale-resources", help="Scale memory/time (env: SCALE_RESOURCES)."
    ),
    verbose: bool = typer.Option(
        False, "--verbose", "-v", help="Increase logging verbosity."
    ),
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
    provided after the known options are passed directly to Snakemake. For example:
        seqnado pipeline rna -n --unlock
    will run a dry-run of the RNA-seq pipeline and unlock the working directory if needed.
    """
    _configure_logging(verbose)

    if show_version:
        from importlib.metadata import version as _pkg_version

        logger.info(f"SeqNado version {_pkg_version('seqnado')}")
        raise typer.Exit(code=0)

    if not _snakemake_available():
        logger.error(
            "`snakemake` not found on PATH. Install/activate the environment that provides it."
        )
        raise typer.Exit(code=127)

    # Local import (keeps single-point-of-truth)
    from seqnado.utils import extract_cores_from_options

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

    # Resolve profile if requested
    profile_trav = None
    if preset:
        profiles = _preset_profiles()
        profile_dir_name = profiles.get(preset.lower())
        if profile_dir_name:
            profile_trav = (
                pkg_root_trav.joinpath("workflow")
                .joinpath("envs")
                .joinpath("profiles")
                .joinpath(profile_dir_name)
            )

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
                if preset and profile_path:
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
                TOP_LEVEL_PASS_THROUGH = (
                    "-n",
                    "--dry-run",
                    "--printshellcmds",
                    "--unlock",
                    "--rerun-incomplete",
                    "--show-failed-logs",
                )

                def should_pass_to_top_level(opt: str) -> bool:
                    # This accepts exact matches and single-token flags that startwith a pass-through key.
                    for p in TOP_LEVEL_PASS_THROUGH:
                        if opt == p or opt.startswith(p + "="):
                            return True
                    # Also accept short flags like "-n" that may be standalone in cleaned_opts list.
                    return False

                # Filter cleaned_opts into top-level opts (preserving order)
                top_level_opts = [o for o in cleaned_opts if should_pass_to_top_level(o)]

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
            if preset and profile_path:
                cmd += ["--profile", str(profile_path)]
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



# -------------------------------- Entrypoint --------------------------------
if __name__ == "__main__":
    app()
