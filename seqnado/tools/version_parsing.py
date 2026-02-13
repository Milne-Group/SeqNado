"""
Version parsing utilities for tool version detection.

This module contains helper functions for extracting version information
from command output, handling apptainer/singularity output formatting,
and parsing version strings.
"""

import re
from typing import List, Optional


def split_output_lines(output: str) -> List[str]:
    """Split output into non-empty lines."""
    return [line for line in output.split("\n") if line.strip()]


def filter_apptainer_info(lines: List[str]) -> List[str]:
    """Remove Apptainer INFO: lines from output."""
    return [line for line in lines if not line.lstrip().startswith("INFO:")]


def first_non_info_line(output: str) -> Optional[str]:
    """Get the first non-INFO line from output."""
    lines = filter_apptainer_info(split_output_lines(output))
    if lines:
        return lines[0]
    raw_lines = split_output_lines(output)
    return raw_lines[0] if raw_lines else None


def extract_version_line(output: str) -> Optional[str]:
    """
    Extract the line containing version information from command output.

    Attempts to intelligently parse output to find version numbers while
    filtering out help text, usage information, and error messages.

    Args:
        output: Command output string

    Returns:
        Line containing version information, or None if not found
    """
    lines = filter_apptainer_info(split_output_lines(output))
    if not lines:
        lines = split_output_lines(output)

    version_pattern = re.compile(
        r"\b(?:v)?\d+(?:\.\d+)+(?:[a-z0-9\-]*)\b", re.IGNORECASE
    )

    for line in lines:
        lower_line = line.lower()

        # Skip lines that are clearly help/usage text
        if any(
            skip in lower_line
            for skip in [
                "usage:",
                "usage example",
                "options:",
                "arguments:",
                "-h",
                "--help",
                "see ",
                "try ",
                "error",
                "invalid",
                "not a valid",
                "unrecognized",
                "command",
                "positional",
                "subcommand",
                "for more information",
            ]
        ):
            continue

        # Skip lines with angle brackets or too many dashes (likely help syntax)
        if "<" in line or ">" in line or line.count("-") > 5:
            continue

        # Look for version keyword followed by a version number
        if "version" in lower_line:
            match = version_pattern.search(line)
            if match:
                return line

        # Or just find a version number pattern on its own line
        match = version_pattern.search(line)
        if match:
            # Make sure it's not embedded in help text
            if not any(
                skip in lower_line
                for skip in ["example", "usage", "type", "default", "options"]
            ):
                return line

    return None


def extract_version_number(output: str) -> Optional[str]:
    """
    Extract just the version number from command output.

    Removes tool names and other text, returning only the version number.

    Args:
        output: Command output string (can be full output or a single line)

    Returns:
        Clean version number string, or None if not found

    Examples:
        >>> extract_version_number("bamnado 0.4.4")
        "0.4.4"
        >>> extract_version_number("FastQC v0.12.1")
        "0.12.1"
        >>> extract_version_number("version 2.5.4")
        "2.5.4"
    """
    # First get the line with version info
    line = extract_version_line(output)
    if not line:
        return None

    # Pattern to match version numbers (with or without 'v' prefix)
    version_pattern = re.compile(
        r"\b(?:v)?(\d+(?:\.\d+)+(?:[a-z0-9\-]*)?)\b", re.IGNORECASE
    )

    match = version_pattern.search(line)
    if match:
        # Return just the version number (group 1), without 'v' prefix
        return match.group(1)

    return None
