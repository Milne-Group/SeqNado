"""SeqNado Command Line Interface package.

This package contains the CLI application and related helper modules.
Commands are being gradually refactored into this package structure.

Public API:
    app: The main Typer application
    
Helper modules (for internal use by commands):
    utils: Logging, file I/O, package resources
    autocomplete: Typer autocomplete functions
    data_validation: Pandera schema validation and data coercion
    snakemake_builder: Centralized Snakemake command construction
    multiomics_helpers: Multiomics mode detection and configuration
    tools_helpers: Tools command helper functions
"""

from __future__ import annotations

# Import the main app from app.py (which contains all the CLI commands)
from seqnado.cli.app import app

# Re-export helper functions for backward compatibility with tests and other code
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
