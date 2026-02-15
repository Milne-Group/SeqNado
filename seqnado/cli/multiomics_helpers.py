"""CLI-specific helpers for multiomics mode detection and configuration."""

from __future__ import annotations

from pathlib import Path
from typing import Optional

from loguru import logger
from seqnado.outputs.multiomics import find_assay_config_paths, find_metadata_paths


def should_use_multiomics_mode(
    assay: Optional[str],
    config_file: Optional[Path],
    cwd: Path = Path("."),
) -> bool:
    """
    Determine if multiomics mode should be used.
    
    Returns True only if:
    - No assay explicitly provided
    - No config file explicitly provided
    - Multiple config_*.yaml files exist in current directory
    
    Args:
        assay: Explicitly provided assay name (if any)
        config_file: Explicitly provided config file path (if any)
        cwd: Working directory to check (default: current directory)
    
    Returns:
        True if multiomics mode should be activated, False otherwise
    """
    # If explicitly provided, don't use multiomics
    if config_file or assay:
        return False

    # Check if multiple config files exist
    config_paths = find_assay_config_paths(cwd)
    return len(config_paths) > 1


def get_multiomics_config_paths(cwd: Path = Path(".")) -> dict:
    """
    Get assay-specific config and metadata paths for multiomics mode.
    
    Args:
        cwd: Working directory to search (default: current directory)
    
    Returns:
        Dictionary with 'configs' and 'metadata' keys, each mapping Assay -> Path
    """
    config_paths = find_assay_config_paths(cwd)
    metadata_paths = find_metadata_paths(cwd)
    return {"configs": config_paths, "metadata": metadata_paths}


def validate_multiomics_setup(cwd: Path = Path(".")) -> bool:
    """
    Validate that all required configs and metadata exist for multiomics mode.
    
    Args:
        cwd: Working directory to validate (default: current directory)
    
    Returns:
        True if setup is valid, False otherwise (logs errors)
    """
    from seqnado.outputs.multiomics import validate_config_and_metadata

    try:
        config_paths = find_assay_config_paths(cwd)
        metadata_paths = find_metadata_paths(cwd)
        validate_config_and_metadata(config_paths, metadata_paths)
        return True
    except FileNotFoundError as e:
        logger.error(str(e))
        return False
