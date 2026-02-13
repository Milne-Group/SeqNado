"""
Container version retrieval strategies for SeqNado tools.

This module handles extracting package version information from various
container types using different strategies (conda-meta, pip list, etc.).
"""

from __future__ import annotations

import json
import re
import subprocess
from typing import Any, Dict, List, Optional

from loguru import logger

# Universal strategy order for extracting package versions from containers
# Ordered from fastest to slowest - tries each until one succeeds
CONTAINER_STRATEGIES: List[Dict[str, Any]] = [
    {
        "type": "conda-meta",
        "description": "Fastest: Check common conda-meta directory (SeqNado containers)",
        "paths": ["/opt/conda/conda-meta"],
    },
    {
        "type": "conda-meta",
        "description": "Fast: Check biocontainer conda-meta directory",
        "paths": ["/usr/local/conda-meta"],
    },
    {
        "type": "pip-list",
        "description": "Medium: Query pip packages (works for Python-based containers)",
    },
    {
        "type": "version-command",
        "description": "Slowest: Run individual tool --version commands (fallback only)",
    },
]


# Mapping from tool name (in tools.json) to conda/pip package name
# for tools whose package name differs from their tool name.
CONDA_PACKAGE_MAP: Dict[str, str] = {
    "bedToBigBed": "ucsc-bedtobigbed",
    "fasterq-dump": "sra-tools",
    "featureCounts": "subread",
    "findPeaks": "homer",
    "macs": "macs2",
    "makeTagDirectory": "homer",
    "ucsc-tools": "ucsc-bedgraphtobigwig",
}

# Module-level cache: container URI -> {package_name: version}
_container_version_cache: Dict[str, Dict[str, str]] = {}

# Conda-meta directories to try (micromamba vs biocontainer layouts)
CONDA_META_DIRS = ["/opt/conda/conda-meta", "/usr/local/conda-meta"]

# Pattern for conda-meta filenames: {name}-{version}-{build}.json
CONDA_META_RE = re.compile(r"^(.+)-([^-]+)-[^-]+\.json$")


def get_conda_package_candidates(
    tool_name: str, tool_info: Optional[Dict[str, Any]] = None
) -> List[str]:
    """Return candidate conda/pip package names for a tool, in priority order.

    Args:
        tool_name: Name of the tool
        tool_info: Optional tool info dict (to avoid circular import)

    Returns:
        List of candidate package names
    """
    if tool_name in CONDA_PACKAGE_MAP:
        return [CONDA_PACKAGE_MAP[tool_name]]

    candidates = [tool_name.lower()]
    if tool_info:
        cmd = tool_info.get("command", "").lower()
        if cmd and cmd not in candidates:
            candidates.append(cmd)
    return candidates


def get_container_strategy(container_uri: str) -> List[Dict[str, Any]]:
    """Get version retrieval strategy for a container.

    Returns the universal ordered strategy list that works from fastest to slowest:
    1. conda-meta at /opt/conda/conda-meta (fastest, common for SeqNado containers)
    2. conda-meta at /usr/local/conda-meta (fast, common for biocontainers)
    3. pip list (medium speed, works for Python-based tools)
    4. version-command (slowest, individual tool queries)

    Args:
        container_uri: Full container URI or image name (currently unused,
                      kept for API compatibility)

    Returns:
        List of strategy configurations to try in order
    """
    return CONTAINER_STRATEGIES


def parse_conda_meta(ls_output: str) -> Dict[str, str]:
    """Parse ``ls`` output of a conda-meta directory into {name: version}."""
    versions: Dict[str, str] = {}
    for filename in ls_output.splitlines():
        m = CONDA_META_RE.match(filename.strip())
        if m:
            versions[m.group(1).lower()] = m.group(2)
    return versions


def parse_pip_json(pip_output: str) -> Dict[str, str]:
    """Parse ``pip list --format=json`` output into {name: version}."""
    packages = json.loads(pip_output)
    return {pkg["name"].lower(): pkg["version"] for pkg in packages}


def get_versions_via_conda_meta(
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
                versions = parse_conda_meta(result.stdout)
                if versions:
                    return versions
        except (subprocess.TimeoutExpired, Exception):
            continue
    return None


def get_versions_via_pip_list(
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
            versions = parse_pip_json(result.stdout)
            if versions:
                return versions
    except (subprocess.TimeoutExpired, json.JSONDecodeError, Exception):
        pass
    return None


def get_versions_via_version_command(
    apptainer_cmd: str, container_uri: str, tool_name: str = None
) -> Optional[Dict[str, str]]:
    """Try to get version from tool's --version command.

    This is useful for containers where metadata is not available.
    Falls back to querying pip or conda if available.
    """
    # If no tool specified, try pip list as fallback
    if tool_name is None:
        return get_versions_via_pip_list(apptainer_cmd, container_uri)

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
            versions = parse_pip_json(result.stdout)
            if versions:
                return versions
    except (subprocess.TimeoutExpired, json.JSONDecodeError, Exception):
        pass

    return None


def get_container_versions(
    container_uri: Optional[str],
    apptainer_cmd: Optional[str],
) -> Dict[str, str]:
    """Get all package versions installed in a container.

    Tries strategies in order from fastest to slowest:

    1. **conda-meta at /opt/conda/conda-meta**: Fastest, used by SeqNado containers
    2. **conda-meta at /usr/local/conda-meta**: Fast, used by biocontainers
    3. **pip list**: Medium speed, works for Python-based containers
    4. **version-command**: Slowest, runs individual tool commands (rarely used)

    The first successful strategy is used and results are cached per container URI.

    Args:
        container_uri: Container URI to query
        apptainer_cmd: Apptainer/singularity command to use

    Returns:
        Mapping of ``{package_name: version}`` (names lower-cased).
    """
    if not container_uri or not apptainer_cmd:
        return {}

    if container_uri in _container_version_cache:
        return _container_version_cache[container_uri]

    # Get strategy for this container
    strategies = get_container_strategy(container_uri)

    # Try each strategy in order
    for strategy in strategies:
        strategy_type = strategy.get("type")

        if strategy_type == "conda-meta":
            paths = strategy.get(
                "paths", ["/opt/conda/conda-meta", "/usr/local/conda-meta"]
            )
            versions = get_versions_via_conda_meta(apptainer_cmd, container_uri, paths)
            if versions:
                _container_version_cache[container_uri] = versions
                logger.debug(
                    f"Retrieved from conda-meta: {container_uri}"
                )
                return versions

        elif strategy_type == "pip-list":
            versions = get_versions_via_pip_list(apptainer_cmd, container_uri)
            if versions:
                _container_version_cache[container_uri] = versions
                logger.debug(
                    f"Retrieved from pip list: {container_uri}"
                )
                return versions

        elif strategy_type == "version-command":
            versions = get_versions_via_version_command(apptainer_cmd, container_uri)
            if versions:
                _container_version_cache[container_uri] = versions
                logger.debug(
                    f"Retrieved from version command: {container_uri}"
                )
                return versions

    logger.debug(f"Retrieved via fallback strategy: {container_uri}")
    _container_version_cache[container_uri] = {}
    return {}


def clear_version_cache():
    """Clear the container version cache. Useful for testing."""
    _container_version_cache.clear()


def extract_version_from_container_tag(container_uri: str) -> Optional[str]:
    """Extract version from container tag/URI.

    Handles various container URI formats:
    - docker://quay.io/biocontainers/meme:5.5.9--pl5321h1ca524f_0 -> "5.5.9"
    - oras://ghcr.io/owner/image:v1.2.3 -> "1.2.3"
    - library://owner/collection/image:1.0 -> "1.0"

    Args:
        container_uri: Full container URI with tag

    Returns:
        Version string extracted from tag, or None if not parseable
    """
    if not container_uri or ":" not in container_uri:
        return None

    # Extract tag after last colon (handles docker://registry:5000/image:tag)
    parts = container_uri.rsplit(":", 1)
    if len(parts) != 2:
        return None

    tag = parts[1]

    # Skip generic tags
    if tag.lower() in ["latest", "stable", "main", "master"]:
        return None

    # Try to extract version-like patterns
    # Pattern 1: Version at start (e.g., "5.5.9--pl5321h1ca524f_0" -> "5.5.9")
    import re

    version_pattern = re.compile(r"^v?(\d+(?:\.\d+)+(?:\.\d+)?)")
    match = version_pattern.match(tag)
    if match:
        return match.group(1)

    # Pattern 2: Just the tag if it looks like a version (e.g., "1.2.3" or "v1.2.3")
    if re.match(r"^v?\d+(?:\.\d+)+$", tag):
        return tag.lstrip("v")

    return None
