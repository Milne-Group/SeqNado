import os
import shutil
import subprocess
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).parent.parent
_DOCKER_IMAGE = "seqnado-devcontainer:test"
# Python for the in-container seqnado venv. Pinned because the image's system
# Python can be newer than some dependencies publish wheels for (e.g.
# `immutables` has none for 3.14, and the image has no Python headers to build it).
_VENV_PYTHON = "3.12"
# The venv volume holds one venv per Python version plus the uv-managed Python
# it links to, so both survive between container runs.
_VENV_ROOT = "/tmp/seqnado_venv"
_VENV_PATH = f"{_VENV_ROOT}/py{_VENV_PYTHON}"
_UV_PYTHON_DIR = f"{_VENV_ROOT}/python"
_APPTAINER_CACHE_PATH = "/apptainer-cache"
# Named Docker volumes persisting the Apptainer image cache and the seqnado
# venv between runs, mapped to their mount points inside the container.
_DOCKER_VOLUMES = {
    "seqnado-apptainer-cache": _APPTAINER_CACHE_PATH,
    "seqnado-test-venv": _VENV_ROOT,
}
_DOCKER_FALLBACK_PRESET = "ld"


def _has_apptainer_runtime() -> bool:
    """Return whether a local Apptainer/Singularity runtime is available."""
    return bool(shutil.which("apptainer") or shutil.which("singularity"))


def _has_docker_runtime() -> bool:
    """Return whether Docker is available for the local_docker test fallback."""
    return bool(shutil.which("docker"))


def pytest_sessionstart(session: pytest.Session) -> None:
    """Ensure a clean pytest basetemp directory for this repo."""
    test_output = Path(__file__).parent.parent / "test_output"
    if test_output.exists():
        shutil.rmtree(test_output, ignore_errors=True)


@pytest.fixture
def sample_metadata_dict():
    """Provide sample metadata dictionary for testing."""
    return {
        "scaling_group": "treatment",
        "condition": "high_dose",
        "spikein_ratio": 0.25,
    }


@pytest.fixture
def multiple_metadata_dicts():
    """Provide multiple metadata dictionaries for testing."""
    return [
        {"scaling_group": "control", "condition": "baseline"},
        {"scaling_group": "treatment", "condition": "low_dose", "spikein_ratio": 0.1},
        {"scaling_group": "treatment", "condition": "high_dose", "spikein_ratio": 0.25},
    ]


@pytest.fixture
def sample_names():
    """Provide sample names for testing."""
    return ["sample1", "sample2", "sample3", "control1", "control2"]


def pytest_configure(config):
    """Configure pytest with custom markers and settings."""
    config.addinivalue_line("markers", "unit: mark test as a unit test")
    config.addinivalue_line("markers", "integration: mark test as an integration test")
    config.addinivalue_line("markers", "edge_case: mark test as an edge case test")
    config.addinivalue_line("markers", "slow: mark test as slow running")
    config.addinivalue_line(
        "markers", "requires_data: mark test as requiring data files"
    )
    config.addinivalue_line(
        "markers", "pipeline: end-to-end Snakemake pipeline tests (opt-in)"
    )
    config.addinivalue_line(
        "markers", "snakemake: tests that invoke Snakemake via subprocess"
    )
    config.addinivalue_line(
        "markers", "genome_build: tests for the genome-build workflow (skipped by default)"
    )
    config.addinivalue_line(
        "markers", "requires_apptainer: tests that require Apptainer/Singularity"
    )


def pytest_addoption(parser: pytest.Parser) -> None:
    """CLI options to control pipeline tests.

    --run-pipeline: Opt-in flag to execute slow pipeline tests.
    --assays: Comma-separated list of assays to run for pipeline tests.
    --cores: Number of cores to use when running the pipeline.
    """
    group = parser.getgroup("seqnado")
    group.addoption(
        "--run-pipeline",
        action="store_true",
        default=False,
        help="Run slow pipeline tests (disabled by default)",
    )
    group.addoption(
        "--assays",
        action="store",
        default="chip",
        help=(
            "Comma-separated list of assays to test "
            "(choices: atac, chip, chip-rx, rna, rna-rx, snp, cat, crispr, meth, mcc). "
            "Default: chip"
        ),
    )
    group.addoption(
        "--cores",
        action="store",
        type=int,
        default=2,
        help="Number of CPU cores to allocate to pipeline tests (default: 2)",
    )
    group.addoption(
        "--run-genome-tests",
        action="store_true",
        default=False,
        help="Run genome-build pipeline tests (disabled by default)",
    )
    group.addoption(
        "--preset",
        action="store",
        default="t",
        help=(
            "Snakemake profile preset for pipeline tests "
            "(e.g. t, tc, lc, ls, le, ld). Default: t"
        ),
    )


def pytest_collection_modifyitems(config: pytest.Config, items: list[pytest.Item]):
    """Modify test collection to add markers based on test names and locations."""
    preset = config.getoption("--preset", default="t")
    apptainer_available = _has_apptainer_runtime()
    docker_available = _has_docker_runtime()

    for item in items:
        # Add unit marker to all tests by default
        if not any(
            mark.name in ["integration", "slow"] for mark in item.iter_markers()
        ):
            item.add_marker(pytest.mark.unit)

        # Add edge_case marker to tests with "edge" in name
        if "edge" in item.name.lower():
            item.add_marker(pytest.mark.edge_case)

        # Add requires_data marker to tests that use test data
        if any(
            fixture in item.fixturenames
            for fixture in ["test_data_dir", "test_fastq_dir"]
        ):
            item.add_marker(pytest.mark.requires_data)

        # If tests are marked as genome_build and user did not opt in, skip them
        if any(m.name == "genome_build" for m in item.iter_markers()):
            if not config.getoption("--run-genome-tests"):
                item.add_marker(
                    pytest.mark.skip(
                        reason="Use --run-genome-tests to enable genome build tests"
                    )
                )

        # If tests are marked as pipeline and user did not opt in, skip them
        elif any(m.name == "pipeline" for m in item.iter_markers()):
            if not config.getoption("--run-pipeline"):
                item.add_marker(
                    pytest.mark.skip(
                        reason="Use --run-pipeline to enable pipeline tests"
                    )
                )

        # Tests that use containerized Snakemake profiles need either a local
        # runtime or the Docker-based local_docker fallback runner.
        if any(m.name == "requires_apptainer" for m in item.iter_markers()):
            uses_docker_fallback = (
                preset == _DOCKER_FALLBACK_PRESET
                and "seqnado_runner" in item.fixturenames
                and docker_available
            )
            if not apptainer_available and not uses_docker_fallback:
                item.add_marker(
                    pytest.mark.skip(
                        reason=(
                            "requires Apptainer/Singularity, or Docker with "
                            "--preset=ld for Docker-backed pipeline tests"
                        )
                    )
                )


# -------------------------
# Global lightweight fixtures
# -------------------------


def pytest_generate_tests(metafunc: pytest.Metafunc) -> None:
    """Parametrize tests that accept the 'assay' fixture using --assays values."""
    if "assay" in metafunc.fixturenames:
        config = metafunc.config
        assays_str: str = config.getoption("--assays") or "chip"
        assays = [a.strip() for a in assays_str.split(",") if a.strip()]
        metafunc.parametrize("assay", assays or ["chip"])


@pytest.fixture(scope="session")
def docker_image(pytestconfig):
    """Build the devcontainer Docker image once per session (for --preset=ld on Mac).

    The image is only built when --preset=ld is specified and Apptainer is not
    available locally (i.e. macOS without Apptainer).
    """
    preset = pytestconfig.getoption("--preset", default="t")
    apptainer_available = _has_apptainer_runtime()
    if preset != "ld" or apptainer_available:
        return None

    if not _has_docker_runtime():
        pytest.skip("Docker not found — required for --preset=ld without Apptainer")

    result = subprocess.run(
        [
            "docker", "build",
            "--platform", "linux/amd64",
            "-t", _DOCKER_IMAGE,
            "-f", str(_REPO_ROOT / ".devcontainer" / "Dockerfile"),
            str(_REPO_ROOT),
        ],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        pytest.fail(f"Failed to build Docker image:\n{result.stderr}")

    # Docker creates new named volumes owned by root, but the image runs as an
    # unprivileged user, so hand the volume mount points over to that user.
    volume_args = [arg for name, path in _DOCKER_VOLUMES.items() for arg in ("-v", f"{name}:{path}")]
    result = subprocess.run(
        [
            "docker", "run", "--rm",
            "--platform", "linux/amd64",
            *volume_args,
            "--entrypoint", "bash",
            "--user", "root",
            _DOCKER_IMAGE,
            "-c",
            f'owner="$(id -u seqnado_user):$(id -g seqnado_user)" '
            f"&& mkdir -p {_APPTAINER_CACHE_PATH}/tmp "
            f'&& chown "$owner" {" ".join(_DOCKER_VOLUMES.values())} {_APPTAINER_CACHE_PATH}/tmp',
        ],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        pytest.fail(f"Failed to prepare Docker volumes:\n{result.stderr}")
    return _DOCKER_IMAGE


def _home_env_args(mounted_root: Path) -> list[str]:
    """Forward HOME into the container when it lives under a mounted path."""
    home = os.environ.get("HOME")
    if home and Path(home).resolve().is_relative_to(mounted_root):
        return ["-e", f"HOME={Path(home).resolve()}"]
    return []


@pytest.fixture
def seqnado_runner(pytestconfig, docker_image, tmp_path_factory):
    """Return a function that runs a seqnado command.

    With --preset=ld and no local Apptainer, each call is routed through
    ``docker run`` using the pre-built devcontainer image.  ``cwd`` and the
    repo root are mounted at their exact host paths so that absolute paths
    embedded in config files remain valid inside the container.
    """
    preset = pytestconfig.getoption("--preset", default="t")
    apptainer_available = _has_apptainer_runtime()
    use_docker = preset == "ld" and not apptainer_available and docker_image is not None

    # Test data (genome, FASTQs) lives in pytest's base temp dir beside each
    # run directory, so mount the whole base temp dir, not just ``cwd``.
    base_temp = tmp_path_factory.getbasetemp().resolve()

    def run(cmd, cwd, **kwargs):
        cwd = Path(cwd).resolve()
        if not use_docker:
            return subprocess.run(cmd, cwd=cwd, **kwargs)

        # Install seqnado into a persistent venv volume (created on first call).
        install = (
            f"export UV_PYTHON_INSTALL_DIR={_UV_PYTHON_DIR} "
            f"&& ({_VENV_PATH}/bin/python -V >/dev/null 2>&1 "
            f"|| uv venv --clear --python {_VENV_PYTHON} {_VENV_PATH}) "
            f"&& . {_VENV_PATH}/bin/activate "
            # Editable, matching how developers run the pipeline.
            f"&& uv pip install -q -e {_REPO_ROOT}"
        )
        # The container has real Apptainer, so run the test profile there: the
        # pipeline helpers add the test-data bind mounts to `profile_test` only.
        cmd = [str(a) for a in cmd]
        for i, arg in enumerate(cmd[:-1]):
            if arg == "--preset" and cmd[i + 1] == _DOCKER_FALLBACK_PRESET:
                cmd[i + 1] = "t"
        seqnado_cmd = " ".join(cmd)

        docker_cmd = [
            "docker", "run", "--rm",
            "--privileged",
            "--platform", "linux/amd64",
            # Mount cwd and repo at identical paths — keeps all absolute
            # references in config files and Apptainer bind-args valid
            "-v", f"{base_temp}:{base_temp}",
            *([] if cwd.is_relative_to(base_temp) else ["-v", f"{cwd}:{cwd}"]),
            "-v", f"{_REPO_ROOT}:{_REPO_ROOT}",
            # Persist the Apptainer image cache and seqnado venv between runs
            *[arg for name, path in _DOCKER_VOLUMES.items() for arg in ("-v", f"{name}:{path}")],
            "-e", f"APPTAINER_CACHEDIR={_APPTAINER_CACHE_PATH}",
            "-e", f"APPTAINER_TMPDIR={_APPTAINER_CACHE_PATH}/tmp",
            "-e", f"SINGULARITY_CACHEDIR={_APPTAINER_CACHE_PATH}",
            "-e", f"SINGULARITY_TMPDIR={_APPTAINER_CACHE_PATH}/tmp",
            # Pipeline fixtures point HOME at the run directory, where `seqnado
            # init` installs the profiles (incl. the bind-patched profile_test).
            *_home_env_args(base_temp),
            "-w", str(cwd),
            docker_image,
            "bash", "-c", f"{install} && {seqnado_cmd}",
        ]
        return subprocess.run(docker_cmd, **kwargs)

    return run
