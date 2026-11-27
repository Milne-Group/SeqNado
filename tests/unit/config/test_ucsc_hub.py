"""Tests for UCSCHubConfig per-assay defaults.

The per-assay grouping defaults must only reference metadata fields that the
hub metadata extractor actually produces for that assay - otherwise tracknado's
TrackDesign crashes with a missing-column assertion at hub build time.
"""

from pathlib import Path

import pytest

from seqnado.core import Assay
from seqnado.config.configs import UCSCHubConfig
from seqnado.workflow.helpers.hub import extract_seqnado_track_metadata

# A representative set of real output paths per assay, covering signal, peaks,
# and the comparison/merged variants so every referenced column is exercised.
REPRESENTATIVE_PATHS = {
    Assay.ATAC: [
        "seqnado_output/atac/bigwigs/deeptools/cpm/sampleA.bigWig",
        "seqnado_output/atac/peaks/lanceotron/sampleA.bed",
    ],
    Assay.CHIP: [
        "seqnado_output/chip/bigwigs/deeptools/cpm/donor-rep1_H3K27ac.bigWig",
        "seqnado_output/chip/peaks/macs2/donor-rep1_H3K27ac.bed",
        "seqnado_output/chip/bigwigs/deeptools/aggregated/treated.bigWig",
    ],
    Assay.CAT: [
        "seqnado_output/cat/bigwigs/bamnado/cpm/s1_CTCF.bigWig",
        "seqnado_output/cat/peaks/seacr/s1_CTCF.bed",
    ],
    Assay.RNA: [
        "seqnado_output/rna/bigwigs/star/cpm/geneA_plus.bigWig",
        "seqnado_output/rna/bigwigs/star/cpm/geneA_minus.bigWig",
    ],
    Assay.MCC: [
        "seqnado_output/mcc/bigwigs/mcc/replicates/sampleA_Myc.bigWig",
        "seqnado_output/mcc/bigwigs/mcc/aggregated-using-mean/groupA_Myc.bigWig",
    ],
    Assay.SNP: [
        "seqnado_output/snp/bigwigs/deeptools/cpm/sampleA.bigWig",
    ],
}


def _referenced_columns(cfg: UCSCHubConfig) -> set[str]:
    cols: set[str] = set()
    for group in (cfg.supergroup_by, cfg.subgroup_by, cfg.overlay_by, cfg.color_by):
        if group:
            cols.update(group)
    return cols


@pytest.mark.parametrize("hub_assay", list(REPRESENTATIVE_PATHS))
def test_defaults_reference_only_producible_fields(hub_assay):
    cfg = UCSCHubConfig.for_assay(hub_assay)
    referenced = _referenced_columns(cfg)

    for path in REPRESENTATIVE_PATHS[hub_assay]:
        produced = set(
            extract_seqnado_track_metadata(Path(path), hub_assay.clean_name)
        )
        missing = referenced - produced
        assert not missing, (
            f"{hub_assay.name} hub default groups on {missing} which the extractor "
            f"does not produce for {path}"
        )


def test_supergroup_and_subgroup_are_disjoint():
    # tracknado asserts supergroup and subgroup columns do not overlap.
    for assay in Assay:
        try:
            cfg = UCSCHubConfig.for_assay(assay)
        except ValueError:
            continue
        supers = set(cfg.supergroup_by or [])
        subs = set(cfg.subgroup_by or [])
        assert supers.isdisjoint(subs), f"{assay.name}: {supers & subs} in both"


def test_rna_overlays_strands_and_colours_by_strand():
    cfg = UCSCHubConfig.for_assay(Assay.RNA)
    assert "strand" in (cfg.color_by or [])
    assert "samplename" in (cfg.overlay_by or [])


def test_chip_sections_by_antibody():
    cfg = UCSCHubConfig.for_assay(Assay.CHIP)
    assert cfg.supergroup_by is not None and "antibody" in cfg.supergroup_by
