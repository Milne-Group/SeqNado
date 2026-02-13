"""
Tools management for SeqNado.

This module provides functionality to list and interact with bioinformatics tools
available in the SeqNado pipeline container.
"""

from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from loguru import logger

from .container_strategies import (
    extract_version_from_container_tag,
    get_conda_package_candidates,
)
from .container_strategies import (
    get_container_versions as _get_container_versions,
)
from .version_parsing import (
    extract_version_line,
    extract_version_number,
    filter_apptainer_info,
    split_output_lines,
)


# Load tools from JSON file at module import time
def _load_tools() -> Dict[str, Dict[str, Any]]:
    """Load tool definitions from JSON file."""
    json_path = Path(__file__).parent.joinpath("tools.json")
    try:
        data = json.loads(json_path.read_text())
        # Extract tools from nested structure
        return data.get("tools", data)
    except Exception as e:
        logger.error(f"Failed to load tools from {json_path}: {e}")
        return {}


AVAILABLE_TOOLS = _load_tools()


def get_container_versions(container_uri: Optional[str] = None) -> Dict[str, str]:
    """Get all package versions installed in a container.

    Tries strategies in order from fastest to slowest:

    1. **conda-meta at /opt/conda/conda-meta**: Fastest, used by SeqNado containers
    2. **conda-meta at /usr/local/conda-meta**: Fast, used by biocontainers
    3. **pip list**: Medium speed, works for Python-based containers
    4. **version-command**: Slowest, runs individual tool commands (rarely used)

    The first successful strategy is used and results are cached per container URI.

    Returns:
        Mapping of ``{package_name: version}`` (names lower-cased).
    """
    if container_uri is None:
        container_uri = get_seqnado_container_uri()
    if not container_uri:
        return {}

    apptainer_cmd = get_apptainer_command()
    if not apptainer_cmd:
        return {}

    return _get_container_versions(container_uri, apptainer_cmd)


def get_tool_version_from_container(tool_name: str) -> Optional[str]:
    """Look up a tool's version from its container's conda/pip metadata.

    Uses :func:`get_container_versions` (cached per container) so that dozens
    of lookups only trigger one ``apptainer exec`` call per unique container.

    Priority order:
    1. Explicit version field in tools.json
    2. Container metadata (conda/pip packages)
    3. Container tag version
    4. None

    Returns:
        Version string, or ``None`` if the package was not found.
    """
    tool_info = get_tool_info(tool_name)
    if not tool_info:
        return None

    # Priority 1: Check for explicit version in tools.json
    if "version" in tool_info and tool_info["version"]:
        return tool_info["version"]

    container_uri = tool_info.get("container") or get_seqnado_container_uri()

    # Skip container detection if tool runs from environment
    if not container_uri or container_uri == "Environment":
        return None

    # Priority 2: Try to get version from container metadata (conda/pip)
    versions = get_container_versions(container_uri)
    if versions:
        for candidate in get_conda_package_candidates(tool_name, tool_info):
            if candidate in versions:
                return versions[candidate]

    # Priority 3: Extract version from container tag
    tag_version = extract_version_from_container_tag(container_uri)
    if tag_version:
        return tag_version

    return None


def _run_local_command(cmd: List[str], timeout: int) -> Tuple[int, str, str]:
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        timeout=timeout,
    )
    return result.returncode, result.stdout, result.stderr


def get_available_tools() -> Dict[str, Dict[str, Any]]:
    """Get dictionary of all available tools."""
    return AVAILABLE_TOOLS


def list_tools(category: Optional[str] = None) -> List[tuple]:
    """
    List available tools, optionally filtered by category.

    Args:
        category: Optional category to filter tools by

    Returns:
        List of tuples (tool_name, description, category)
    """
    # Get categories in preferred order
    all_categories = get_categories()
    category_order = {cat: i for i, cat in enumerate(all_categories)}

    tools = []

    for tool_name, tool_info in AVAILABLE_TOOLS.items():
        cat = tool_info["category"]
        # Handle both single string and list of categories
        cats = cat if isinstance(cat, list) else [cat]

        if category is None or category.lower() in [c.lower() for c in cats]:
            # For display, use first category or all if multiple
            display_cat = cats[0] if cats else "Unknown"
            tools.append((tool_name, tool_info["description"], display_cat))

    # Sort by category order, then by tool name
    return sorted(tools, key=lambda x: (category_order.get(x[2], 999), x[0]))


def get_categories() -> List[str]:
    """Get list of all available tool categories in workflow order."""
    # Define preferred category order
    category_order = [
        "Download",
        "Quality Control",
        "Preprocessing",
        "Alignment",
        "Analysis",
        "Visualization",
        "Reporting",
        "Quantification",
        "Utilities",
    ]

    # Collect all unique categories from tools
    categories = set()
    for tool_info in AVAILABLE_TOOLS.values():
        cat = tool_info["category"]
        cats = cat if isinstance(cat, list) else [cat]
        for c in cats:
            categories.add(c)

    # Sort by preferred order, then alphabetically for any missing
    sorted_cats = []
    for cat in category_order:
        if cat in categories:
            sorted_cats.append(cat)
            categories.remove(cat)

    # Add any remaining categories not in the preferred order
    sorted_cats.extend(sorted(categories))

    return sorted_cats


def tool_exists(tool_name: str) -> bool:
    """Check if a tool exists in the available tools list."""
    return tool_name.lower() in {name.lower() for name in AVAILABLE_TOOLS.keys()}


def get_tool_info(tool_name: str) -> Optional[Dict[str, Any]]:
    """Get information about a specific tool."""
    for name, info in AVAILABLE_TOOLS.items():
        if name.lower() == tool_name.lower():
            return info
    return None


def _get_tool_subcommands(tool_info: Dict[str, Any]) -> List[str]:
    subcommands = tool_info.get("subcommands", [])
    if isinstance(subcommands, str):
        subcommands = [subcommands]
    if not isinstance(subcommands, list):
        return []
    return [str(cmd) for cmd in subcommands if str(cmd).strip()]


def get_tool_subcommands(tool_name: str) -> List[str]:
    """Get subcommands for a tool, if defined."""
    tool_info = get_tool_info(tool_name)
    if not tool_info:
        return []
    return _get_tool_subcommands(tool_info)


def _resolve_tool_command(tool_info: Dict[str, Any], subcommand: Optional[str]) -> str:
    if not subcommand:
        return tool_info["command"]

    subcommands = _get_tool_subcommands(tool_info)
    for cmd in subcommands:
        if cmd.lower() == subcommand.lower():
            return cmd

    return tool_info["command"]


def get_tool_version(
    tool_name: str,
    use_container: bool = False,
    subcommand: Optional[str] = None,
) -> str:
    """
    Get version information for a tool.

    Attempts multiple strategies in order:
    1. If using container, query container metadata first (cleanest version)
    2. Run the tool with common version flags (``--version``, ``-v``)
       either locally or inside the container
    3. Fall back to "Version information not available"

    Args:
        tool_name: Name of the tool
        use_container: Whether to run the command in a container
        subcommand: Optional subcommand

    Returns:
        Version string or error message
    """
    tool_info = get_tool_info(tool_name)
    if not tool_info:
        return f"Tool '{tool_name}' not found"

    # Priority 1: Try container metadata first for clean version numbers
    if use_container:
        version = get_tool_version_from_container(tool_name)
        if version:
            return version

    # Check if tool runs from environment rather than container
    container_uri = tool_info.get("container")
    runs_from_environment = not container_uri or container_uri == "Environment"

    # Priority 2: Try running --version command
    command = _resolve_tool_command(tool_info, subcommand)

    version_flags = [
        ["--version"],
        ["-v"],
    ]

    for flags in version_flags:
        try:
            # Use container only if use_container=True AND tool doesn't run from environment
            if use_container and not runs_from_environment:
                returncode, stdout, stderr = run_command_in_container(
                    tool_name,
                    flags,
                    subcommand=subcommand,
                )
                output = (stdout + stderr).strip()
                if output and returncode in [0, 1]:
                    version = extract_version_number(output)
                    if version:
                        return version
                continue

            cmd = [command] + flags
            returncode, stdout, stderr = _run_local_command(cmd, timeout=5)
            output = (stdout + stderr).strip()
            if output and returncode in [0, 1]:
                version = extract_version_number(output)
                if version:
                    return version
        except (FileNotFoundError, subprocess.TimeoutExpired, Exception):
            continue

    return "Version information not available"


def get_tool_help(
    tool_name: str,
    use_container: bool = False,
    subcommand: Optional[str] = None,
) -> str:
    """
    Get help information for a tool.

    Args:
        tool_name: Name of the tool
        use_container: Whether to run the command in a container

    Returns:
        Help output as string
    """
    tool_info = get_tool_info(tool_name)
    if not tool_info:
        return f"Tool '{tool_name}' not found"

    command = _resolve_tool_command(tool_info, subcommand)

    # Check if tool runs from environment rather than container
    container_uri = tool_info.get("container")
    runs_from_environment = not container_uri or container_uri == "Environment"

    try:
        # Try to get help with common help flags
        help_flags = ["--help", "-h", "-help", "help"]

        for flag in help_flags:
            try:
                # Use container only if use_container=True AND tool doesn't run from environment
                if use_container and not runs_from_environment:
                    returncode, stdout, stderr = run_command_in_container(
                        tool_name,
                        [flag],
                        subcommand=subcommand,
                    )
                    output = (stdout + stderr).strip()
                    if output and returncode in [0, 1]:
                        lines = split_output_lines(output)
                        filtered = filter_apptainer_info(lines)
                        return "\n".join(filtered if filtered else lines)
                    continue

                cmd = [command, flag]
                returncode, stdout, stderr = _run_local_command(cmd, timeout=10)
                if returncode == 0 or stdout or stderr:
                    output = stdout or stderr
                    if output.strip():
                        return output.strip()
            except Exception:
                continue

        return f"Help information not available for '{command}'"

    except Exception as e:
        logger.debug(f"Error getting help for {tool_name}: {e}")
        return f"Error retrieving help for '{tool_name}': {str(e)}"


def _get_bibtex_entry(citation_key: str) -> Optional[str]:
    """Extract a raw BibTeX entry from the citations file by key."""
    bib_path = Path(__file__).parent.joinpath("tool_citations.bib")
    if not bib_path.exists():
        return None

    try:
        content = bib_path.read_text()
    except Exception as e:
        logger.debug(f"Failed to read citations file: {e}")
        return None

    pattern = re.compile(
        rf"^(@\w+\{{{re.escape(citation_key)},.*?^}})",
        re.MULTILINE | re.DOTALL,
    )
    match = pattern.search(content)
    if match:
        return match.group(1)

    return None


def format_citation(bibtex: str) -> Optional[str]:
    """
    Format a BibTeX entry into a human-readable citation string.

    Args:
        bibtex: Raw BibTeX entry string

    Returns:
        Formatted citation string, or None if formatting fails
    """
    try:
        import io

        from pybtex.backends.plaintext import Backend
        from pybtex.database import parse_string
        from pybtex.style.formatting.plain import Style

        bib_data = parse_string(bibtex, "bibtex")
        style = Style()
        formatted = style.format_bibliography(bib_data)

        backend = Backend()
        output = io.StringIO()
        for entry in formatted:
            output.write(entry.text.render(backend))

        result = output.getvalue().strip()
        return result if result else None
    except Exception as e:
        logger.debug(f"Failed to format citation: {e}")
        return None


def get_tool_citation(tool_name: str) -> Optional[str]:
    """
    Get the BibTeX citation for a tool, if available.

    Args:
        tool_name: Name of the tool

    Returns:
        BibTeX entry as a string, or None if no citation is available
    """
    tool_info = get_tool_info(tool_name)
    if not tool_info:
        return None

    citation_key = tool_info.get("citation")
    if not citation_key:
        return None

    return _get_bibtex_entry(citation_key)


def check_tool_available(tool_name: str, subcommand: Optional[str] = None) -> bool:
    """
    Check if a tool is available on the system.

    Args:
        tool_name: Name of the tool

    Returns:
        True if tool is available, False otherwise
    """
    tool_info = get_tool_info(tool_name)
    if not tool_info:
        return False

    command = _resolve_tool_command(tool_info, subcommand)
    return shutil.which(command) is not None


def get_seqnado_container_uri() -> Optional[str]:
    """
    Extract the container URI from the SeqNado Snakefile.

    Returns:
        Container URI string (e.g., "oras://ghcr.io/alsmith151/seqnado_pipeline:latest")
        or None if not found
    """
    try:
        from importlib import resources

        snakefile_trav = (
            resources.files("seqnado").joinpath("workflow").joinpath("Snakefile")
        )

        with resources.as_file(snakefile_trav) as snakefile_path:
            content = Path(snakefile_path).read_text()

            # Look for: container: "URI"
            for line in content.split("\n"):
                line = line.strip()
                if line.startswith("container:"):
                    # Extract the URI from quotes
                    # Format: container: "oras://ghcr.io/..."
                    if '"' in line:
                        uri = line.split('"')[1]
                        return uri
    except Exception as e:
        logger.debug(f"Could not extract container URI: {e}")

    return None


def is_apptainer_available() -> bool:
    """Check if apptainer is available on the system."""
    return (
        shutil.which("apptainer") is not None or shutil.which("singularity") is not None
    )


def get_apptainer_command() -> Optional[str]:
    """Get the apptainer/singularity command to use."""
    if shutil.which("apptainer"):
        return "apptainer"
    elif shutil.which("singularity"):
        return "singularity"
    return None


def run_command_in_container(
    tool_name: str,
    args: List[str],
    container_uri: Optional[str] = None,
    subcommand: Optional[str] = None,
) -> Tuple[int, str, str]:
    """
    Run a command in the appropriate container for the tool.

    Args:
        tool_name: Name of the tool
        args: Arguments to pass to the tool
        container_uri: Container URI (defaults to tool-specific container from tools.json)

    Returns:
        Tuple of (return_code, stdout, stderr)
    """
    tool_info = get_tool_info(tool_name)
    if not tool_info:
        return 1, "", f"Tool '{tool_name}' not found"

    # Use provided container URI, or get from tool info, or fall back to default
    if container_uri is None:
        container_uri = tool_info.get("container")

    if container_uri is None:
        container_uri = get_seqnado_container_uri()

    if not container_uri:
        return 1, "", "Could not determine container URI"

    apptainer_cmd = get_apptainer_command()
    if not apptainer_cmd:
        return 1, "", "Apptainer/Singularity not available on system"

    command = _resolve_tool_command(tool_info, subcommand)

    env = os.environ.copy()

    # Build the command to run in container
    cmd = [apptainer_cmd, "exec", container_uri, command] + args

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=30,
            env=env,
        )
        return result.returncode, result.stdout, result.stderr
    except subprocess.TimeoutExpired:
        return 1, "", "Command timed out after 30 seconds"
    except Exception as e:
        return 1, "", f"Error running command: {str(e)}"


def run_tool_help_in_container(
    tool_name: str,
    container_uri: Optional[str] = None,
    subcommand: Optional[str] = None,
) -> str:
    """
    Get help for a tool by running it in the container.

    Args:
        tool_name: Name of the tool
        container_uri: Container URI (defaults to seqnado pipeline container)

    Returns:
        Help output as string
    """
    tool_info = get_tool_info(tool_name)
    if not tool_info:
        return f"Tool '{tool_name}' not found"

    # Try common help flags
    help_flags_list = [
        ["--help"],
        ["-h"],
    ]

    for help_flags in help_flags_list:
        returncode, stdout, stderr = run_command_in_container(
            tool_name,
            help_flags,
            container_uri=container_uri,
            subcommand=subcommand,
        )

        # Combine stdout and stderr
        output = (stdout + stderr).strip()
        if output and returncode in [0, 1]:  # Some tools return 1 for help
            lines = split_output_lines(output)
            filtered = filter_apptainer_info(lines)
            return "\n".join(filtered if filtered else lines)

    subcommands = _get_tool_subcommands(tool_info)
    if subcommands and not subcommand:
        return (
            "This tool exposes subcommands. Use --subcommand to request help for one "
            f"of: {', '.join(subcommands)}"
        )

    return f"Help information not available for '{_resolve_tool_command(tool_info, subcommand)}'"
