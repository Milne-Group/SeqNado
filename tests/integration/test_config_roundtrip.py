from __future__ import annotations

import os
import subprocess
from pathlib import Path

import pytest

from seqnado.config.core import SeqnadoConfig
from tests.pipeline.helpers.data import GenomeResources
from tests.pipeline.helpers import data as pipeline_data


def _mock_rna_resources(genome_path: Path) -> GenomeResources:
    genome_path.mkdir(parents=True, exist_ok=True)

    bt2_dir = genome_path / "bt2_chr21"
    bt2_dir.mkdir()
    for suffix in GenomeResources._BT2_SUFFIXES:
        (bt2_dir / f"bt2_chr21.{suffix}").write_text("index\n")

    star_dir = genome_path / "STAR_chr21_rna_spikein"
    star_dir.mkdir()

    files = {
        "chromosome_sizes": genome_path / "chr21.chrom.sizes",
        "gtf": genome_path / "chr21_rna_spikein.gtf",
        "blacklist": genome_path / "hg38_chr21-blacklist.bed",
        "genes": genome_path / "hg38_genes.bed",
        "fasta": genome_path / "chr21.fa",
        "plot_coords": genome_path / "plotting_coordinates.bed",
    }
    for path in files.values():
        path.write_text("test\n")

    return GenomeResources(
        assay="rna",
        star_index=star_dir,
        bt2_index=bt2_dir / "bt2_chr21",
        chromosome_sizes=files["chromosome_sizes"],
        gtf=files["gtf"],
        blacklist=files["blacklist"],
        genes=files["genes"],
        fasta=files["fasta"],
        plot_coords=files["plot_coords"],
    )


@pytest.mark.integration
def test_cli_config_roundtrip_with_tests_data(tmp_path: Path):
    genome_path = tmp_path / "genome"

    assay = "rna"
    resources = _mock_rna_resources(genome_path)

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


@pytest.mark.integration
def test_download_resources_reuses_existing_bt2_index_when_refresh_fails(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    genome_path = tmp_path / "genome"
    index_dir = genome_path / "bt2_chr21"
    index_dir.mkdir(parents=True)

    for suffix in GenomeResources._BT2_SUFFIXES:
        (index_dir / f"bt2_chr21.{suffix}").write_text("index\n")

    real_has_complete = GenomeResources._has_complete_bt2_index
    calls = {"count": 0}

    def flaky_has_complete(index_dir: Path, prefix: str) -> bool:
        calls["count"] += 1
        if calls["count"] == 1:
            return False
        return real_has_complete(index_dir, prefix)

    monkeypatch.setattr(GenomeResources, "_has_complete_bt2_index", flaky_has_complete)
    monkeypatch.setattr(
        pipeline_data,
        "download_with_retry",
        lambda *args, **kwargs: (_ for _ in ()).throw(FileNotFoundError("missing tar")),
    )

    result = GenomeResources._ensure_index(
        genome_path, "https://example.invalid", "bt2_chr21", prefix="bt2_chr21"
    )

    assert result == index_dir / "bt2_chr21"
