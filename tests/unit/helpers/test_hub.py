"""Tests for seqnado.workflow.helpers.hub metadata extraction."""

from pathlib import Path
from types import SimpleNamespace

import pytest
import pandas as pd

from seqnado.workflow.helpers.hub import (
    METADATA_KEYS,
    extract_seqnado_track_metadata,
    format_hub_labels,
    load_design_metadata,
    make_seqnado_extractor,
)

PREFIX = "seqnado_output"


def meta(path: str, assay: str) -> dict[str, str]:
    return extract_seqnado_track_metadata(Path(f"{PREFIX}/{path}"), assay)


class TestSchema:
    @pytest.mark.parametrize(
        "track_path,assay_name",
        [
            ("bigwigs/deeptools/cpm/sampleA.bigWig", "atac"),
            ("peaks/macs2/sampleA_H3K4me3.bed", "chip"),
            ("bigwigs/star/cpm/geneA_plus.bigWig", "rna"),
            ("bigwigs/mcc/replicates/sampleA_vpX.bigWig", "mcc"),
        ],
    )
    def test_all_keys_present(self, track_path, assay_name):
        m = meta(track_path, assay_name)
        assert set(m) == set(METADATA_KEYS)
        # No key is ever missing; unused fields are "NA".
        assert all(v is not None for v in m.values())


class TestAtac:
    def test_signal(self):
        m = meta("bigwigs/deeptools/cpm/sampleA.bigWig", "atac")
        assert m["file_type"] == "signal"
        assert m["method"] == "deeptools"
        assert m["norm"] == "cpm"
        assert m["samplename"] == "sampleA"
        assert m["track_kind"] == "individual"
        assert m["antibody"] == "NA"

    def test_peak(self):
        m = meta("peaks/lanceotron/sampleA.bed", "atac")
        assert m["file_type"] == "peaks"
        assert m["method"] == "lanceotron"
        assert m["norm"] == "NA"  # peaks have no scaling level
        assert m["samplename"] == "sampleA"


class TestChip:
    def test_antibody_from_signal(self):
        m = meta("bigwigs/deeptools/cpm/donor-rep1_H3K27ac.bigWig", "chip")
        assert m["samplename"] == "donor-rep1"
        assert m["antibody"] == "H3K27ac"
        assert m["norm"] == "cpm"

    def test_antibody_from_peak(self):
        m = meta("peaks/macs2/donor-rep1_H3K27ac.bed", "chip")
        assert m["file_type"] == "peaks"
        assert m["antibody"] == "H3K27ac"
        assert m["samplename"] == "donor-rep1"

    def test_spikein(self):
        m = meta("bigwigs/deeptools/spikein/orlando/donor-rep1_H3K27ac.bigWig", "chip")
        assert m["norm"] == "spikein"
        assert m["spikein"] == "orlando"
        assert m["method"] == "deeptools"
        assert m["antibody"] == "H3K27ac"

    def test_cat_alias(self):
        m = meta("bigwigs/bamnado/cpm/s1_CTCF.bigWig", "cat")
        assert m["antibody"] == "CTCF"

    def test_comparison_has_no_antibody(self):
        # Aggregated/subtraction tracks are condition-level, not per-antibody.
        agg = meta("bigwigs/deeptools/aggregated/treated.bigWig", "chip")
        assert agg["track_kind"] == "aggregated"
        assert agg["antibody"] == "NA"
        assert agg["samplename"] == "treated"

        sub = meta("bigwigs/deeptools/subtraction/treated_vs_control.bigWig", "chip")
        assert sub["track_kind"] == "subtraction"
        assert sub["antibody"] == "NA"
        assert sub["samplename"] == "treated_vs_control"


class TestRna:
    @pytest.mark.parametrize("strand", ["plus", "minus"])
    def test_strand(self, strand):
        m = meta(f"bigwigs/star/cpm/geneA_{strand}.bigWig", "rna")
        assert m["strand"] == strand
        assert m["samplename"] == "geneA"

    def test_no_antibody_parsing_for_rna(self):
        m = meta("bigwigs/star/cpm/geneA_plus.bigWig", "rna")
        assert m["antibody"] == "NA"


class TestMcc:
    def test_replicate(self):
        m = meta("bigwigs/mcc/replicates/sampleA_Myc.bigWig", "mcc")
        assert m["method"] == "mcc"
        assert m["viewpoint"] == "Myc"
        assert m["samplename"] == "sampleA"
        assert m["track_kind"] == "individual"

    def test_aggregated_keeps_viewpoint(self):
        m = meta("bigwigs/mcc/aggregated-using-mean/groupA_Myc.bigWig", "mcc")
        assert m["track_kind"] == "aggregated"
        assert m["viewpoint"] == "Myc"

    def test_subtraction(self):
        m = meta("bigwigs/mcc/subtractions/g1_vs_g2_Myc.bigWig", "mcc")
        assert m["track_kind"] == "subtraction"
        assert m["viewpoint"] == "Myc"
        assert m["samplename"] == "g1_vs_g2"


class TestMisc:
    def test_merged(self):
        m = meta("bigwigs/deeptools/merged/csaw/groupA.bigWig", "chip")
        assert m["merged"] == "yes"
        assert m["norm"] == "csaw"

    def test_bigbed_suffix(self):
        m = meta("peaks/macs2/sampleA_H3K4me3.bb", "chip")
        assert m["file_type"] == "peaks"
        assert m["antibody"] == "H3K4me3"

    def test_unknown_layout_still_builds(self):
        m = extract_seqnado_track_metadata(Path("weird/place/track.bigWig"), "atac")
        assert m["samplename"] == "track"
        assert set(m) == set(METADATA_KEYS)

    def test_make_extractor_binds_assay(self):
        extractor = make_seqnado_extractor("rna")
        m = extractor(Path(f"{PREFIX}/bigwigs/star/cpm/geneA_minus.bigWig"))
        assert m["assay"] == "rna"
        assert m["strand"] == "minus"


class TestDesignMetadata:
    def test_loads_standard_and_user_defined_columns(self, tmp_path):
        design = tmp_path / "metadata.csv"
        design.write_text(
            "sample_id,uid,condition,donor\n"
            "donor1_rep1,donor1_rep1_H3K27ac,treated,D1\n"
        )

        design_metadata = load_design_metadata(design)
        m = extract_seqnado_track_metadata(
            Path(f"{PREFIX}/bigwigs/deeptools/cpm/donor1_rep1_H3K27ac.bigWig"),
            "chip",
            design_metadata,
        )

        assert m["condition"] == "treated"
        assert m["donor"] == "D1"
        assert m["antibody"] == "H3K27ac"

    def test_design_metadata_fields_are_present_for_unmatched_tracks(self, tmp_path):
        design = tmp_path / "metadata.csv"
        design.write_text("sample_id,condition\nsampleA,treated\n")

        m = extract_seqnado_track_metadata(
            Path(f"{PREFIX}/bigwigs/deeptools/aggregated/treated.bigWig"),
            "chip",
            load_design_metadata(design),
        )

        assert m["condition"] == "NA"


class TestHubLabels:
    def test_supertrack_labels_are_concise_and_descriptive(self):
        composite = SimpleNamespace(
            tracktype="bigWig",
            kwargs={},
        )
        supertrack = SimpleNamespace(
            subtracks=[composite],
            kwargs={},
        )
        hub = SimpleNamespace(
            track_design=SimpleNamespace(
                details=pd.DataFrame(
                    [{
                        "supertrack": "group-1",
                        "file_type": "signal",
                        "condition": "treated",
                        "antibody": "H3K27ac",
                    }]
                ),
                super_tracks={"group-1": supertrack},
            )
        )

        format_hub_labels(hub, ["file_type", "condition", "antibody"])

        assert supertrack.kwargs == {
            "shortLabel": "treated H3K27ac",
            "longLabel": "File Type: signal · Condition: treated · Antibody: H3K27ac",
        }
        assert composite.kwargs == {
            "shortLabel": "Signal tracks",
            "longLabel": "Signal tracks — File Type: signal · Condition: treated · Antibody: H3K27ac",
        }
