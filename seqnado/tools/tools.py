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


def _split_output_lines(output: str) -> List[str]:
    return [line for line in output.split("\n") if line.strip()]


def _filter_apptainer_info(lines: List[str]) -> List[str]:
    return [line for line in lines if not line.lstrip().startswith("INFO:")]


def _first_non_info_line(output: str) -> Optional[str]:
    lines = _filter_apptainer_info(_split_output_lines(output))
    if lines:
        return lines[0]
    raw_lines = _split_output_lines(output)
    return raw_lines[0] if raw_lines else None


def _extract_version_line(output: str) -> Optional[str]:
    lines = _filter_apptainer_info(_split_output_lines(output))
    if not lines:
        lines = _split_output_lines(output)

    version_pattern = re.compile(r"\b(?:v)?\d+(?:\.\d+)+(?:[a-z0-9\-]*)\b", re.IGNORECASE)
    
    for line in lines:
        lower_line = line.lower()
        
        # Skip lines that are clearly help/usage text
        if any(skip in lower_line for skip in [
            "usage:", "usage example", "options:", "arguments:", "-h", "--help",
            "see ", "try ", "error", "invalid", "not a valid", "unrecognized",
            "command", "positional", "subcommand", "for more information"
        ]):
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
            if not any(skip in lower_line for skip in ["example", "usage", "type", "default", "options"]):
                return line

    return None


# Mapping from tool name (in tools.json) to conda/pip package name
# for tools whose package name differs from their tool name.
_CONDA_PACKAGE_MAP: Dict[str, str] = {
    "bedToBigBed": "ucsc-bedtobigbed",
    "fasterq-dump": "sra-tools",
    "featureCounts": "subread",
    "findPeaks": "homer",
    "makeTagDirectory": "homer",
    "macs": "macs2",
    "lanceotron-mcc": "lanceotron",
    "ucsc-tools": "ucsc-bedgraphtobigwig",
}

# Module-level cache: container URI -> {package_name: version}
_container_version_cache: Dict[str, Dict[str, str]] = {}


def _get_conda_package_candidates(tool_name: str) -> List[str]:
    """Return candidate conda/pip package names for a tool, in priority order."""
    if tool_name in _CONDA_PACKAGE_MAP:
        return [_CONDA_PACKAGE_MAP[tool_name]]

    candidates = [tool_name.lower()]
    tool_info = get_tool_info(tool_name)
    if tool_info:
        cmd = tool_info.get("command", "").lower()
        if cmd and cmd not in candidates:
            candidates.append(cmd)
    return candidates


# Conda-meta directories to try (micromamba vs biocontainer layouts)
_CONDA_META_DIRS = ["/opt/conda/conda-meta", "/usr/local/conda-meta"]

# Pattern for conda-meta filenames: {name}-{version}-{build}.json
_CONDA_META_RE = re.compile(r"^(.+)-([^-]+)-[^-]+\.json$")


# Container-specific version retrieval strategies
# Maps container identifier to strategy configuration
_CONTAINER_STRATEGIES: Dict[str, Dict[str, Any]] = {
    # SeqNado pipeline (micromamba conda)
    "ghcr.io/alsmith151/seqnado_pipeline": {
        "strategies": [
            {"type": "conda-meta", "paths": ["/opt/conda/conda-meta"]},
            {"type": "pip-list"},
        ]
    },
    # BAM processing container (needs --version)
    "ghcr.io/alsmith151/bamnado": {
        "strategies": [
            {"type": "version-command"},
            {"type": "pip-list"},
        ]
    },
    # MCC analysis container (pip-only)
    "ghcr.io/alsmith151/mccnado": {
        "strategies": [
            {"type": "pip-list"},
            {"type": "conda-meta", "paths": ["/opt/conda/conda-meta", "/usr/local/conda-meta"]},
        ]
    },
    # Biocontainers (conda in /usr/local/conda-meta)
    "quay.io/biocontainers": {
        "strategies": [
            {"type": "conda-meta", "paths": ["/usr/local/conda-meta"]},
            {"type": "pip-list"},
        ]
    },
    # SeqNado ML CPU (pip-only)
    "ghcr.io/alsmith151/seqnado_ml_cpu": {
        "strategies": [
            {"type": "pip-list"},
        ]
    },
    # PlotNado (unknown - try standard fallbacks)
    "asmith151/plotnado": {
        "strategies": [
            {"type": "conda-meta", "paths": ["/opt/conda/conda-meta", "/usr/local/conda-meta"]},
            {"type": "pip-list"},
        ]
    },
}


def _get_container_strategy(container_uri: str) -> List[Dict[str, Any]]:
    """Get version retrieval strategy for a container.
    
    Args:
        container_uri: Full container URI or image name
    
    Returns:
        List of strategy configurations to try in order
    """
    # Try to match against known patterns
    for pattern, config in _CONTAINER_STRATEGIES.items():
        if pattern in container_uri:
            return config.get("strategies", [])
    
    # Default fallback strategy
    return [
        {"type": "conda-meta", "paths": ["/opt/conda/conda-meta", "/usr/local/conda-meta"]},
        {"type": "pip-list"},
    ]


def _parse_conda_meta(ls_output: str) -> Dict[str, str]:
    """Parse ``ls`` output of a conda-meta directory into {name: version}."""
    versions: Dict[str, str] = {}
    for filename in ls_output.splitlines():
        m = _CONDA_META_RE.match(filename.strip())
        if m:
            versions[m.group(1).lower()] = m.group(2)
    return versions


def _parse_pip_json(pip_output: str) -> Dict[str, str]:
    """Parse ``pip list --format=json`` output into {name: version}."""
    packages = json.loads(pip_output)
    return {pkg["name"].lower(): pkg["version"] for pkg in packages}


def _get_versions_via_conda_meta(
    apptainer_cmd: str, container_uri: str, paths: List[str]
) -> Optional[Dict[str, str]]:
    """Try to get versions from conda-meta directory."""
    for meta_dir in paths:
        try:
            result = subprocess.run(
                [apptainer_cmd, "exec", container_uri, "ls", meta_dir],
                capture_output=True,
                text=True,
                timeout=30,
            )
            if result.returncode == 0 and result.stdout.strip():
                versions = _parse_conda_meta(result.stdout)
                if versions:
                    return versions
        except (subprocess.TimeoutExpired, Exception):
            continue
    return None


def _get_versions_via_pip_list(
    apptainer_cmd: str, container_uri: str
) -> Optional[Dict[str, str]]:
    """Try to get versions from pip list."""
    try:
        result = subprocess.run(
            [apptainer_cmd, "exec", container_uri, "pip", "list", "--format=json"],
            capture_output=True,
            text=True,
            timeout=30,
        )
        if result.returncode == 0 and result.stdout.strip():
            versions = _parse_pip_json(result.stdout)
            if versions:
                return versions
    except (subprocess.TimeoutExpired, json.JSONDecodeError, Exception):
        pass
    return None


def _get_versions_via_version_command(
    apptainer_cmd: str, container_uri: str, tool_name: str = None
) -> Optional[Dict[str, str]]:
    """Try to get version from tool's --version command.
    
    This is useful for containers where metadata is not available.
    Falls back to querying pip or conda if available.
    """
    # If no tool specified, try pip list as fallback
    if tool_name is None:
        return _get_versions_via_pip_list(apptainer_cmd, container_uri)
    
    try:
        # List available tools in the container by looking at what's in PATH
        # or try common package managers
        result = subprocess.run(
            [apptainer_cmd, "exec", container_uri, "pip", "list", "--format=json"],
            capture_output=True,
            text=True,
            timeout=30,
        )
        if result.returncode == 0 and result.stdout.strip():
            versions = _parse_pip_json(result.stdout)
            if versions:
                return versions
    except (subprocess.TimeoutExpired, json.JSONDecodeError, Exception):
        pass
    
    return None


def get_container_versions(container_uri: Optional[str] = None) -> Dict[str, str]:
    """Get all package versions installed in a container.

    Attempts container-specific strategies based on known container types:
    
    - **seqnado_pipeline**: conda-meta at /opt/conda/conda-meta, then pip
    - **bamnado**: pip list (no conda metadata available)
    - **mccnado**: pip list, then conda-meta fallback
    - **biocontainers** (mageck, meme): conda-meta at /usr/local/conda-meta, then pip
    - **seqnado_ml_cpu**: pip list
    - **other/unknown**: Try conda-meta dirs, then pip list

    Results are cached per container URI for the lifetime of the process.

    Returns:
        Mapping of ``{package_name: version}`` (names lower-cased).
    """
    if container_uri is None:
        container_uri = get_seqnado_container_uri()
    if not container_uri:
        return {}

    if container_uri in _container_version_cache:
        return _container_version_cache[container_uri]

    apptainer_cmd = get_apptainer_command()
    if not apptainer_cmd:
        return {}

    # Get strategy for this container
    strategies = _get_container_strategy(container_uri)
    
    # Try each strategy in order
    for strategy in strategies:
        strategy_type = strategy.get("type")
        
        if strategy_type == "conda-meta":
            paths = strategy.get("paths", ["/opt/conda/conda-meta", "/usr/local/conda-meta"])
            versions = _get_versions_via_conda_meta(apptainer_cmd, container_uri, paths)
            if versions:
                _container_version_cache[container_uri] = versions
                logger.debug(f"Retrieved {len(versions)} packages from conda-meta: {container_uri}")
                return versions
        
        elif strategy_type == "pip-list":
            versions = _get_versions_via_pip_list(apptainer_cmd, container_uri)
            if versions:
                _container_version_cache[container_uri] = versions
                logger.debug(f"Retrieved {len(versions)} packages from pip list: {container_uri}")
                return versions
        
        elif strategy_type == "version-command":
            versions = _get_versions_via_version_command(apptainer_cmd, container_uri)
            if versions:
                _container_version_cache[container_uri] = versions
                logger.debug(f"Retrieved {len(versions)} packages from version command: {container_uri}")
                return versions

    logger.debug(f"No package metadata found in {container_uri}")
    _container_version_cache[container_uri] = {}
    return {}


def get_tool_version_from_container(tool_name: str) -> Optional[str]:
    """Look up a tool's version from its container's conda/pip metadata.

    Uses :func:`get_container_versions` (cached per container) so that dozens
    of lookups only trigger one ``apptainer exec`` call per unique container.

    Returns:
        Version string, or ``None`` if the package was not found.
    """
    tool_info = get_tool_info(tool_name)
    if not tool_info:
        return None

    container_uri = tool_info.get("container") or get_seqnado_container_uri()
    if not container_uri:
        return None

    versions = get_container_versions(container_uri)
    if not versions:
        return None

    for candidate in _get_conda_package_candidates(tool_name):
        if candidate in versions:
            return versions[candidate]

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

    Runs the tool with common version flags (``--version``, ``-v``)
    either locally or inside the container.

    For bulk version lookups (e.g. doc generation) prefer
    :func:`get_tool_version_from_container` which queries conda
    metadata once per container image.

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

    command = _resolve_tool_command(tool_info, subcommand)

    version_flags = [
        ["--version"],
        ["-v"],
    ]

    for flags in version_flags:
        try:
            if use_container:
                returncode, stdout, stderr = run_command_in_container(
                    tool_name,
                    flags,
                    subcommand=subcommand,
                )
                output = (stdout + stderr).strip()
                if output and returncode in [0, 1]:
                    line = _extract_version_line(output)
                    if line:
                        return line[:200]
                continue

            cmd = [command] + flags
            returncode, stdout, stderr = _run_local_command(cmd, timeout=5)
            output = (stdout + stderr).strip()
            if output and returncode in [0, 1]:
                line = _extract_version_line(output)
                if line:
                    return line[:200]
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

    try:
        # Try to get help with common help flags
        help_flags = ["--help", "-h", "-help", "help"]

        for flag in help_flags:
            try:
                if use_container:
                    returncode, stdout, stderr = run_command_in_container(
                        tool_name,
                        [flag],
                        subcommand=subcommand,
                    )
                    output = (stdout + stderr).strip()
                    if output and returncode in [0, 1]:
                        lines = _split_output_lines(output)
                        filtered = _filter_apptainer_info(lines)
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
            lines = _split_output_lines(output)
            filtered = _filter_apptainer_info(lines)
            return "\n".join(filtered if filtered else lines)

    subcommands = _get_tool_subcommands(tool_info)
    if subcommands and not subcommand:
        return (
            "This tool exposes subcommands. Use --subcommand to request help for one "
            f"of: {', '.join(subcommands)}"
        )

    return f"Help information not available for '{_resolve_tool_command(tool_info, subcommand)}'"
