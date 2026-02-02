from __future__ import annotations

import os
import subprocess
from pathlib import Path

import pytest

from seqnado.config.core import SeqnadoConfig
from tests.pipeline.helpers.data import GenomeResources


@pytest.mark.integration
def test_cli_config_roundtrip_with_tests_data(tmp_path: Path):
    genome_path = tmp_path / "genome"

    # Download resources for RNA (has STAR index)
    assay = "rna"
    resources = GenomeResources.download_resources(genome_path, assay)

    # Write genome config under isolated HOME/SEQNADO_CONFIG
    os.environ["HOME"] = str(tmp_path)
    os.environ["SEQNADO_CONFIG"] = str(tmp_path)
    genome_config_file = tmp_path / ".config" / "seqnado" / "genome_config.json"
    resources.write_config(genome_config_file)

    # Create seqnado_output dir for UCSCHubConfig validation
    (tmp_path / "seqnado_output").mkdir()

    # Generate config via CLI (non-interactive) into a known file
    out = tmp_path / f"config_{assay}.yaml"
    res = subprocess.run(
        [
            "seqnado",
            "config",
            assay,
            "--no-interactive",
            "--no-make-dirs",
            "--render-options",
            "-o",
            str(out),
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )

    if res.returncode != 0:
        print("STDOUT:\n", res.stdout)
        print("STDERR:\n", res.stderr)
    assert res.returncode == 0
    assert out.exists()

    # Load config into model
    cfg = SeqnadoConfig.from_yaml(out)
    assert cfg.assay.value.lower() == assay
    # Roundtrip: dump to YAML-equivalent dict and reload
    import yaml

    dump = tmp_path / "roundtrip.yaml"
    with dump.open("w") as f:
        yaml.safe_dump(cfg.model_dump(mode="json"), f)
    cfg2 = SeqnadoConfig.from_yaml(dump)
    assert cfg2.assay == cfg.assay
    assert cfg2.project.name == cfg.project.name
