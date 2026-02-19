import shutil
import subprocess
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).parent.parent
_DOCKER_IMAGE = "seqnado-devcontainer:test"
_VENV_PATH = "/tmp/seqnado_venv"

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

        # If tests are marked as pipeline and user did not opt in, skip them
        if any(m.name == "pipeline" for m in item.iter_markers()):
            if not config.getoption("--run-pipeline"):
                item.add_marker(
                    pytest.mark.skip(
                        reason="Use --run-pipeline to enable pipeline tests"
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
    apptainer_available = bool(shutil.which("apptainer") or shutil.which("singularity"))
    if preset != "ld" or apptainer_available:
        return None

    if not shutil.which("docker"):
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
    return _DOCKER_IMAGE


@pytest.fixture
def seqnado_runner(pytestconfig, docker_image):
    """Return a function that runs a seqnado command.

    With --preset=ld and no local Apptainer, each call is routed through
    ``docker run`` using the pre-built devcontainer image.  ``cwd`` and the
    repo root are mounted at their exact host paths so that absolute paths
    embedded in config files remain valid inside the container.
    """
    preset = pytestconfig.getoption("--preset", default="t")
    apptainer_available = bool(shutil.which("apptainer") or shutil.which("singularity"))
    use_docker = preset == "ld" and not apptainer_available and docker_image is not None

    def run(cmd, cwd, **kwargs):
        cwd = Path(cwd).resolve()
        if not use_docker:
            return subprocess.run(cmd, cwd=cwd, **kwargs)

        # Install seqnado into a persistent venv volume on first call;
        # subsequent calls skip the install (uv detects nothing changed).
        install = (
            f"(test -f {_VENV_PATH}/bin/activate || uv venv {_VENV_PATH}) "
            f"&& . {_VENV_PATH}/bin/activate "
            f"&& uv pip install -q -e {_REPO_ROOT}"
        )
        seqnado_cmd = " ".join(str(a) for a in cmd)

        docker_cmd = [
            "docker", "run", "--rm",
            "--privileged",
            "--platform", "linux/amd64",
            # Mount cwd and repo at identical paths — keeps all absolute
            # references in config files and Apptainer bind-args valid
            "-v", f"{cwd}:{cwd}",
            "-v", f"{_REPO_ROOT}:{_REPO_ROOT}",
            # Persist the Apptainer image cache and seqnado venv between runs
            "-v", "seqnado-apptainer-cache:/apptainer-cache",
            "-v", f"seqnado-test-venv:{_VENV_PATH}",
            "-e", "APPTAINER_CACHEDIR=/apptainer-cache",
            "-e", "APPTAINER_TMPDIR=/apptainer-cache/tmp",
            "-e", "SINGULARITY_CACHEDIR=/apptainer-cache",
            "-e", "SINGULARITY_TMPDIR=/apptainer-cache/tmp",
            "-w", str(cwd),
            docker_image,
            "bash", "-c", f"{install} && {seqnado_cmd}",
        ]
        return subprocess.run(docker_cmd, **kwargs)

    return run
