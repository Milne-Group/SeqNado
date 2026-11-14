"""
Unit tests ensuring every container image SeqNado runs is pinned and defined
once in ``seqnado/containers.py``.

Floating tags (``:latest``) are cached by Snakemake under the image URI, so a
newly published image is never re-pulled and the workflow runs against stale
tools. See ``seqnado/containers.py`` for details.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

from seqnado import containers
from seqnado.tools.container_strategies import extract_version_from_container_tag

pytestmark = pytest.mark.unit

PKG_DIR = Path(__file__).parent.parent.parent / "seqnado"
WORKFLOW_DIR = PKG_DIR / "workflow"
URI_RE = re.compile(r"(?:oras|docker|library)://[^\s\"',)\\]+")

WORKFLOW_FILES = sorted(
    [
        *WORKFLOW_DIR.joinpath("rules").rglob("*.smk"),
        *WORKFLOW_DIR.glob("Snakefile*"),
        *WORKFLOW_DIR.joinpath("scripts").rglob("*.py"),
    ]
)
FLOATING_TAGS = {"latest", "stable", "main", "master"}


def _is_pinned(uri: str) -> bool:
    ref = uri.split("://", 1)[1]
    if "@sha256:" in ref or ":sha256." in ref:
        return True
    name = ref.rsplit("/", 1)[-1]
    if ":" not in name:
        return False  # no tag means implicit latest
    return name.rsplit(":", 1)[1].lower() not in FLOATING_TAGS


@pytest.mark.parametrize("name,uri", sorted(containers.ALL_CONTAINERS.items()))
def test_registry_images_are_pinned(name: str, uri: str):
    if uri in containers.UNPINNED_ALLOWLIST:
        pytest.xfail(f"{name} is not yet pinned (see TODO in seqnado/containers.py)")
    assert _is_pinned(uri), f"{name} uses a floating tag: {uri}"


@pytest.mark.parametrize(
    "path", WORKFLOW_FILES, ids=lambda p: str(p.relative_to(WORKFLOW_DIR))
)
def test_workflow_files_do_not_hardcode_images(path: Path):
    """Rules/scripts must reference CONTAINERS.<NAME>, not literal URIs."""
    hardcoded = URI_RE.findall(path.read_text())
    assert not hardcoded, (
        f"{path.name} hardcodes container URI(s) {hardcoded}; "
        "add them to seqnado/containers.py and reference the constant instead"
    )


def test_tools_json_images_are_registered():
    tools = json.loads((PKG_DIR / "tools" / "tools.json").read_text())
    registered = set(containers.ALL_CONTAINERS.values())
    used = {
        info["container"]
        for info in tools["tools"].values()
        if "://" in str(info.get("container", ""))
    }
    assert used, "no container URIs found in tools.json; has its layout changed?"
    assert used <= registered, (
        f"tools.json uses unregistered images: {sorted(used - registered)}"
    )


@pytest.mark.parametrize(
    "uri,expected",
    [
        ("docker://ghcr.io/alsmith151/bamnado:0.9.0", "0.9.0"),
        ("oras://ghcr.io/alsmith151/seqnado_pipeline@sha256:" + "a" * 64, None),
        ("docker://ghcr.io/org/tool:1.2.3@sha256:" + "b" * 64, "1.2.3"),
        ("docker://registry:5000/image@sha256:" + "c" * 64, None),
    ],
)
def test_extract_version_handles_digests(uri: str, expected: str | None):
    assert extract_version_from_container_tag(uri) == expected
