"""SeqNado Command Line Interface package.

This package contains the CLI application and related helper modules.

Public API:
    app: The main Typer application (with all commands registered)
"""

from __future__ import annotations

# Import from app.py (not app_instance.py) to trigger command registration
from seqnado.cli.app import app
from seqnado.cli.autocomplete import (
    _assay_names,
    _find_fastqs,
    _get_profile_name,
    _preset_profiles,
    _profile_autocomplete,
)
from seqnado.cli.data_validation import (
    _coerce_value_to_dtype,
    _extract_candidate_defaults_from_schema,
    _format_col_hint,
)
from seqnado.cli.utils import (
    _configure_logging,
    _pkg_traversable,
    _read_json,
    _snakemake_available,
    _style_name_with_rich,
    _write_json,
)

__all__ = [
    "app",
    "_assay_names",
    "_coerce_value_to_dtype",
    "_configure_logging",
    "_extract_candidate_defaults_from_schema",
    "_find_fastqs",
    "_format_col_hint",
    "_get_profile_name",
    "_pkg_traversable",
    "_preset_profiles",
    "_profile_autocomplete",
    "_read_json",
    "_snakemake_available",
    "_style_name_with_rich",
    "_write_json",
]
