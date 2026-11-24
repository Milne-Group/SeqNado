import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from seqnado import containers as CONTAINERS


@pytest.mark.pipeline
@pytest.mark.snakemake
@pytest.mark.parametrize(
    "container",
    [
        pytest.param(False, id="native"),
        pytest.param(True, id="container", marks=pytest.mark.requires_apptainer),
    ],
)
def test_mcc_peak_rule(tmp_path: Path, container: bool):
    if not container and not shutil.which("mcc-peaks"):
        pytest.skip("mcc-peaks is not installed in the execution environment")

    image = CONTAINERS.MCC_PEAKS
    output_dir = tmp_path / "output files"
    bigwig = output_dir / "bigwigs/mcc/group_vp.bigWig"
    bigwig.parent.mkdir(parents=True)
    viewpoints = tmp_path / "viewpoints.bed"
    viewpoints.write_text("chr1\t0\t10\tvp-chr1-0-10\nchr1\t450\t460\tother\n")
    exclusions = output_dir / "resources/exclusion_regions.bed"
    exclusions.parent.mkdir(parents=True)
    exclusions.write_text("chr1\t100\t120\n")

    if container:
        runtime = shutil.which("apptainer") or shutil.which("singularity")
        assert runtime is not None
        python = [runtime, "exec", image, "python"]
    else:
        python = [sys.executable]
    subprocess.run(
        python
        + [
            "-c",
            (
                "import pybigtools, sys; "
                "w = pybigtools.open(sys.argv[1], 'w'); "
                "w.write({'chr1': 500}, [('chr1', 0, 100, 1.0), "
                "('chr1', 100, 120, 10.0), ('chr1', 120, 300, 1.0), "
                "('chr1', 300, 320, 8.0), ('chr1', 320, 500, 1.0)]); "
                "w.close()"
            ),
            str(bigwig),
        ],
        check=True,
        cwd=tmp_path,
    )

    rule_file = Path(__file__).resolve().parents[2] / "seqnado/workflow/rules/mcc/peaks.smk"
    snakefile = tmp_path / "Snakefile"
    snakefile.write_text(
        "from types import SimpleNamespace as NS\n"
        "from seqnado import containers as CONTAINERS\n"
        "from seqnado.config.third_party_tools import LanceotronMCC\n"
        f"OUTPUT_DIR = {str(output_dir)!r}\n"
        f"viewpoints = {str(viewpoints)!r}\n"
        "tool = LanceotronMCC()\n"
        "tool.command_line_arguments.value += ' --thresholds 1 --background-windows 51"
        " --window-sizes 5 --support-fraction 1 --min-peak-size 3'\n"
        "CONFIG = NS(assay_config=NS(mcc=NS(viewpoints=viewpoints)), "
        "third_party_tools=NS(lanceotron_mcc=tool))\n"
        "SAMPLE_GROUPINGS = NS(get_grouping=lambda name: NS(group_names=['group']))\n"
        "VIEWPOINT_TO_GROUPED_VIEWPOINT = {'vp-chr1-0-10': 'vp'}\n"
        "SCALE_RESOURCES = 1\n"
        "def define_memory_requested(**kwargs): return '1GB'\n"
        "def define_time_requested(**kwargs): return 1\n"
        f"include: {str(rule_file)!r}\n"
    )
    peaks = output_dir / "peaks/lanceotron-mcc/group_vp.bed"
    command = ["snakemake", "--snakefile", str(snakefile), "--cores", "2"]
    if container:
        command.append("--use-apptainer")
    subprocess.run(command + [str(peaks)], cwd=tmp_path, check=True)

    calls = [line.split("\t") for line in peaks.read_text().splitlines()]
    assert any(chrom == "chr1" and int(start) < 320 and int(end) > 300 for chrom, start, end in calls)
    assert all(not (int(start) < 120 and int(end) > 100) for _, start, end in calls)
    assert Path(f"{peaks}.tsv").exists()
    assert Path(f"{peaks}.unresolved.bed").exists()
    metadata = json.loads(Path(f"{peaks}.json").read_text())
    assert metadata["parameters"]["viewpoints"] == {"chr1": [[0, 10]]}
