import json
import re
import tarfile
import time
from dataclasses import dataclass
from pathlib import Path

import pytest
import requests
import yaml


def download_with_retry(url: str, dest: Path, max_retries: int = 3, timeout: int = 30):
    """Download a file with retry logic."""
    for attempt in range(max_retries):
        try:
            dest.parent.mkdir(parents=True, exist_ok=True)
            r = requests.get(url, stream=True, timeout=timeout)
            r.raise_for_status()
            with open(dest, "wb") as f:
                f.write(r.content)
            return
        except (requests.RequestException, OSError):
            if attempt < max_retries - 1:
                time.sleep(2**attempt)
                continue
            raise


def extract_tar(tar_path: Path, dest: Path, flatten: bool = True):
    """Extract a tar archive, optionally flattening nested directories."""
    with tarfile.open(tar_path) as tar:
        tar.extractall(path=dest)

    if flatten:
        nested = dest / dest.name
        if nested.exists() and nested.is_dir():
            for f in nested.iterdir():
                f.rename(dest / f.name)
            try:
                nested.rmdir()
            except Exception as e:
                print(f"[WARNING] Could not remove nested directory {nested}: {e}")


def get_fastq_pattern(assay: str) -> str:
    """Get the FASTQ file pattern for a given assay."""
    patterns = {
        "atac": "atac_*.fastq.gz",
        "cat": "chip-rx_*.fastq.gz",
        "chip-rx": "chip-rx_*.fastq.gz",
        "chip": "chip-rx_*.fastq.gz",
        "crispr": "crispr-*.fastq.gz",
        "mcc": "mcc-*.fastq.gz",
        "meth": "meth-*.fastq.gz",
        "rna-rx": "rna-spikein-*.fastq.gz",
        "rna": "rna*.fastq.gz",
        "snp": "snp_*.fastq.gz",
    }
    pattern = patterns.get(assay)
    if not pattern:
        raise ValueError(f"Unsupported assay: {assay}")
    return pattern


@dataclass(frozen=True)
class TestPaths:
    """Container for all important test and project paths. All paths are created and ready for use."""

    repo: Path
    package: Path
    test_dir: Path
    test_data: Path
    workflow: Path
    test_profile: Path
    genome: Path
    fastq: Path


def make_test_paths(file: Path) -> TestPaths:
    """Factory function for TestPaths."""
    # Find project root by searching upwards for pyproject.toml
    repo = file.resolve()
    while not (repo / "pyproject.toml").exists():
        if repo.parent == repo:
            raise RuntimeError(
                f"Could not find project root (pyproject.toml) from {file}"
            )
        repo = repo.parent

    # Define paths
    package = repo / "seqnado"
    test_dir = repo / "test_output"
    test_data = test_dir / "data"
    genome = test_data / "genome"
    fastq = test_data / "fastq"

    # Create directories
    for path in (test_dir, test_data, genome, fastq):
        path.mkdir(parents=True, exist_ok=True)

    return TestPaths(
        repo=repo,
        package=package,
        test_dir=test_dir,
        test_data=test_data,
        workflow=package / "workflow",
        test_profile=package / "workflow" / "envs" / "profiles" / "profile_test",
        genome=genome,
        fastq=fastq,
    )


class TestContext:
    """Unified context class for test session/function-scoped helpers."""

    def __init__(self, pytestconfig, tmp_path_factory):
        self.pytestconfig = pytestconfig
        self.tmp_path_factory = tmp_path_factory
        self.test_paths = make_test_paths(Path(__file__).resolve())

    @property
    def cores(self):
        return int(self.pytestconfig.getoption("--cores"))

    @property
    def selected_assays(self):
        assays_str: str = self.pytestconfig.getoption("--assays") or "chip"
        assays = [a.strip() for a in assays_str.split(",") if a.strip()]
        return assays or ["chip"]

    def assay_type(self, assay: str) -> str:
        """Extract base assay type from assay name."""
        return re.sub(r"(.*)\-.*", r"\1", assay)

    def run_directory(self, assay: str) -> Path:
        try:
            base_temp = self.tmp_path_factory.getbasetemp()
        except (FileExistsError, AttributeError):
            base_temp = self.tmp_path_factory._basetemp

        # If base_temp is still None, use mktemp to create it
        if base_temp is None:
            base_temp = self.tmp_path_factory.mktemp("pytest")

        run_dir = base_temp / assay
        run_dir.mkdir(exist_ok=True, parents=True)
        return run_dir

    def plot_bed(self, test_data_path: Path):
        return test_data_path / "genome" / "plotting_coordinates.bed"
