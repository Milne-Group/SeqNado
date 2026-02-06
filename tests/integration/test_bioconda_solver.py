from __future__ import annotations

import os
import shutil
import subprocess
import uuid

import pytest


@pytest.mark.integration
@pytest.mark.slow
def test_bioconda_solver_for_seqnado_dependency_set(tmp_path):
    """Regression test for bioconda dependency conflicts reported in #331."""
    if os.getenv("SEQNADO_RUN_BIOCONDA_SOLVER_TEST") != "1":
        pytest.skip("Set SEQNADO_RUN_BIOCONDA_SOLVER_TEST=1 to run bioconda solve test")

    solver = shutil.which("mamba") or shutil.which("micromamba") or shutil.which("conda")
    if not solver:
        pytest.skip("No conda solver found (mamba, micromamba, or conda)")

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
        "numpy",
        "pandas",
        "snakemake-executor-plugin-slurm",
        "tracknado",
        "cookiecutter",
        "pillow",
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
