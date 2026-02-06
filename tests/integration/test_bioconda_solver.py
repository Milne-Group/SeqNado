from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
import uuid
import zipfile
from pathlib import Path

import pytest


def _build_wheel_from_source(output_dir: Path) -> Path:
    """Build a wheel from the current source tree and return its path."""
    output_dir.mkdir(parents=True, exist_ok=True)
    build_cmd = [
        sys.executable,
        "-m",
        "pip",
        "wheel",
        ".",
        "--no-deps",
        "--no-build-isolation",
        "--wheel-dir",
        str(output_dir),
    ]
    env = os.environ.copy()
    env["PIP_CACHE_DIR"] = str(output_dir / ".pip-cache")
    result = subprocess.run(build_cmd, capture_output=True, text=True, env=env)
    assert result.returncode == 0, (
        "Failed to build seqnado wheel from source.\n"
        f"Command: {' '.join(build_cmd)}\n"
        f"STDOUT:\n{result.stdout}\n"
        f"STDERR:\n{result.stderr}"
    )
    wheels = sorted(output_dir.glob("*.whl"))
    assert wheels, f"No wheel produced in {output_dir}"
    return wheels[-1]


def _normalise_requirement_to_conda_spec(requirement: str) -> str:
    """Convert a PEP 508 requirement string to a conda-compatible spec."""
    req = requirement.strip()
    if " @" in req:
        req = req.split(" @", maxsplit=1)[0].strip()

    if " (" in req and req.endswith(")"):
        name, spec = req.split(" (", maxsplit=1)
        spec = spec[:-1]
    else:
        match = re.search(r"[<>=!~]", req)
        if match:
            idx = match.start()
            name, spec = req[:idx], req[idx:]
        else:
            name, spec = req, ""

    name = name.strip().split("[", maxsplit=1)[0]
    spec = spec.replace(" ", "")
    return f"{name}{spec}"


def _metadata_requirements_for_conda(
    wheel_path: Path, include_extras: set[str]
) -> tuple[list[str], str | None]:
    """Extract conda specs and python requirement from wheel metadata."""
    with zipfile.ZipFile(wheel_path, "r") as zip_file:
        metadata_member = next(
            name for name in zip_file.namelist() if name.endswith(".dist-info/METADATA")
        )
        metadata = zip_file.read(metadata_member).decode("utf-8")

    python_spec = None
    conda_specs: list[str] = []
    for line in metadata.splitlines():
        if line.startswith("Requires-Python: "):
            python_spec = line.split(":", maxsplit=1)[1].strip()
            continue

        if not line.startswith("Requires-Dist: "):
            continue

        requirement = line.split(":", maxsplit=1)[1].strip()
        req_part, sep, marker = requirement.partition(";")
        marker = marker.strip().lower()
        if sep and "extra ==" in marker:
            marker = marker.replace("'", '"')
            include = any(
                f'extra == "{extra.lower()}"' in marker for extra in include_extras
            )
            if not include:
                continue

        conda_specs.append(_normalise_requirement_to_conda_spec(req_part))

    return sorted(set(conda_specs)), python_spec


@pytest.mark.integration
@pytest.mark.slow
def test_bioconda_solver_for_seqnado_dependency_set(tmp_path):
    """Build seqnado and verify conda can solve its runtime dependencies."""
    if os.getenv("SEQNADO_RUN_BIOCONDA_SOLVER_TEST") != "1":
        pytest.skip("Set SEQNADO_RUN_BIOCONDA_SOLVER_TEST=1 to run bioconda solve test")

    solver = shutil.which("mamba") or shutil.which("micromamba") or shutil.which("conda")
    if not solver:
        pytest.skip("No conda solver found (mamba, micromamba, or conda)")

    wheel = _build_wheel_from_source(tmp_path / "dist")
    conda_specs, _ = _metadata_requirements_for_conda(wheel, include_extras={"slurm"})
    assert conda_specs, "No runtime requirements extracted from built seqnado wheel"

    env_prefix = tmp_path / f"seqnado-bioconda-solve-{uuid.uuid4().hex[:8]}"
    cmd = [
        solver,
        "create",
        "--yes",
        "--dry-run",
        "--override-channels",
        "-c",
        "conda-forge",
        "-c",
        "bioconda",
        "-p",
        str(env_prefix),
        "python=3.11",
        *conda_specs,
    ]

    env = os.environ.copy()
    # Keep cache writes in a writable location (important on read-only base envs).
    pkgs_dir = tmp_path / "pkgs"
    env["MAMBA_PKGS_DIRS"] = str(pkgs_dir)
    env["CONDA_PKGS_DIRS"] = str(pkgs_dir)

    result = subprocess.run(cmd, capture_output=True, text=True, env=env)
    stderr_lower = result.stderr.lower()
    if (
        "could not resolve host" in stderr_lower
        or "could not resolve hostname" in stderr_lower
        or "temporary failure in name resolution" in stderr_lower
    ):
        pytest.skip("Network unavailable for bioconda solve test")

    assert result.returncode == 0, (
        "Bioconda dependency solve failed for SeqNado conflict set.\n"
        f"Command: {' '.join(cmd)}\n"
        f"STDOUT:\n{result.stdout}\n"
        f"STDERR:\n{result.stderr}"
    )
