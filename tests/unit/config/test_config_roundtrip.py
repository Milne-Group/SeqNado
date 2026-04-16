from __future__ import annotations

from datetime import date
from importlib import resources
from pathlib import Path

import yaml

from seqnado import Assay
from seqnado.config.configs import GenomeConfig, ProjectConfig, STARIndex
from seqnado.config.core import SeqnadoConfig
from seqnado.config.user_input import build_default_assay_config, render_config


def _build_rna_config(tmp_path: Path) -> SeqnadoConfig:
    star_dir = tmp_path / "star"
    star_dir.mkdir()
    (tmp_path / "metadata.csv").write_text("sample_id,fastq\n")

    genome = GenomeConfig(name="hg38", index=STARIndex(prefix=star_dir))
    assay_config = build_default_assay_config(Assay.RNA, genome)
    project = ProjectConfig(name="test_project", date=date.today())

    return SeqnadoConfig(
        assay=Assay.RNA,
        project=project,
        genome=genome,
        metadata=tmp_path / "metadata.csv",
        assay_config=assay_config,
    )


def test_create_bigwigs_roundtrip_toggle(tmp_path: Path) -> None:
    config = _build_rna_config(tmp_path)
    out_file = tmp_path / "config_rna.yaml"
    seqnado_version = "9.9.9"

    template = resources.files("seqnado.data").joinpath("config_template.jinja")
    with resources.as_file(template) as tpl_path:
        render_config(Path(tpl_path), config, out_file, seqnado_version=seqnado_version)

    config_data = yaml.safe_load(out_file.read_text())
    assert "assay_config" in config_data
    assert config_data["seqnado_version"] == seqnado_version

    for value in (True, False):
        config_data["assay_config"]["create_bigwigs"] = value
        out_file.write_text(yaml.safe_dump(config_data, sort_keys=False))
        loaded = SeqnadoConfig.from_yaml(out_file)
        assert loaded.assay_config.create_bigwigs is value
