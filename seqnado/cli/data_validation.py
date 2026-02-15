"""Data validation and schema helpers for design file generation."""

from __future__ import annotations

import re
from typing import TYPE_CHECKING, Any, List, Optional

import typer
from loguru import logger

from seqnado.cli.utils import _style_name_with_rich

if TYPE_CHECKING:
    import pandas as pd


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
) -> tuple[pd.Series, pd.Series] | None:
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
        Tuple of (group_names, binary_encoding) or None if groups cannot be reliably extracted
    """
    import pandas as pd

    # Strategy 1: Regex pattern matching
    if pattern:
        try:
            groups = sample_ids.str.extract(pattern, expand=False)
            if groups.notna().all() and groups.nunique() >= 2:
                # Validate groups are reasonably balanced
                group_counts = groups.value_counts()
                if group_counts.min() > 0 and group_counts.max() / len(groups) <= 0.9:
                    return (groups, None)
        except Exception:
            pass  # Fall through to other strategies

    # Strategy 2 & 3: Keyword detection and heuristics
    # Common group keywords to look for
    group_keywords = [
        "control",
        "ctrl",
        "treated",
        "treatment",
        "treat",
        "wt",
        "wildtype",
        "wild-type",
        "ko",
        "knockout",
        "knock-out",
        "mut",
        "mutant",
        "untreated",
        "knockdown",
        "kd",
        "vehicle",
        "mock",
        "dmso",
    ]

    def extract_group(sample_id: str) -> str:
        """Extract group from a single sample_id."""
        sample_lower = sample_id.lower()

        # Split on both hyphens and underscores
        parts = re.split(r"[-_]", sample_id)

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
    control_keywords = [
        "control",
        "ctrl",
        "untreated",
        "vehicle",
        "mock",
        "dmso",
        "wt",
        "wildtype",
    ]

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
        if (
            accept_all_defaults
            and col in ("deseq2", "group")
            and "sample_id" in df.columns
            and assay is not None
        ):
            from seqnado.inputs import Assay as AssayEnum

            if assay == AssayEnum.RNA:
                result = _extract_deseq2_groups_from_sample_names(
                    df["sample_id"], pattern=deseq2_pattern
                )
                if result is not None:
                    groups, deseq2_binary = result
                    df["group"] = groups
                    if deseq2_binary is not None:
                        df["deseq2"] = deseq2_binary
                        logger.info(
                            f"Adding 'group' and 'deseq2' columns with groups: {sorted(set(groups))}"
                        )
                    else:
                        # 3+ groups case - populate group column and add empty deseq2
                        df["deseq2"] = pd.NA  # Add empty deseq2 column to prevent prompting
                        logger.info(
                            f"Adding 'group' column with {len(set(groups))} groups: {sorted(set(groups))}"
                        )
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
        if (
            col in ("group", "deseq2")
            and "sample_id" in df.columns
            and assay is not None
        ):
            from seqnado.inputs import Assay as AssayEnum

            if assay == AssayEnum.RNA:
                extracted_result = _extract_deseq2_groups_from_sample_names(
                    df["sample_id"], pattern=deseq2_pattern
                )
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
                            logger.info(
                                f"Added 'group' and 'deseq2' columns with auto-detected groups: {unique_groups}"
                            )
                            continue
                    else:
                        # 3+ groups - offer uniform value or use detected groups
                        typer.echo(
                            f"⚠️ Auto-detected {len(unique_groups)} groups: {unique_groups}"
                        )
                        typer.echo(
                            "Multi-group comparisons require manual DESeq2 configuration."
                        )
                        use_detected = typer.confirm(
                            "Use auto-detected groups for 'group' column (leave 'deseq2' empty for manual config)?",
                            default=True,
                        )
                        if use_detected:
                            df["group"] = groups
                            df["deseq2"] = (
                                pd.NA
                            )  # Add empty deseq2 column to prevent prompting
                            logger.info(
                                f"Added 'group' column with {len(unique_groups)} groups"
                            )
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
