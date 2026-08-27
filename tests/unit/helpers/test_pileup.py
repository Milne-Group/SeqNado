"""Tests for seqnado.workflow.helpers.pileup."""

import pytest

from seqnado.core import DataScalingTechnique
from seqnado.inputs.grouping import SampleGroup, SampleGroupings, SampleGroups
from seqnado.workflow.helpers.pileup import get_condition_input_bigwigs


class _Wildcards:
    def __init__(self, condition: str):
        self.condition = condition


@pytest.fixture
def groupings() -> SampleGroupings:
    return SampleGroupings(
        groupings={
            "condition": SampleGroups(
                groups=[SampleGroup(name="treated", samples=["s1", "s2"])]
            )
        }
    )


@pytest.mark.unit
class TestGetConditionInputBigwigs:
    """The condition-aggregation rules feed these paths straight into bamnado, so a
    template that fails to interpolate takes the whole DAG down."""

    def test_per_sample_paths_are_fully_interpolated(self, groupings):
        paths = get_condition_input_bigwigs(
            _Wildcards("treated"),
            pileup_method="deeptools",
            output_dir="out",
            sample_groupings=groupings,
        )
        scale = DataScalingTechnique.PER_SAMPLE.value
        assert paths == [
            f"out/bigwigs/deeptools/{scale}/s1.bigWig",
            f"out/bigwigs/deeptools/{scale}/s2.bigWig",
        ]

    def test_stranded_paths_carry_the_strand_suffix(self, groupings):
        """RNA-seq is the only caller that passes a strand, and the suffix has to
        survive interpolation rather than being left as a literal placeholder."""
        paths = get_condition_input_bigwigs(
            _Wildcards("treated"),
            pileup_method="deeptools",
            output_dir="out",
            sample_groupings=groupings,
            strand="plus",
        )
        scale = DataScalingTechnique.PER_SAMPLE.value
        assert paths == [
            f"out/bigwigs/deeptools/{scale}/s1_plus.bigWig",
            f"out/bigwigs/deeptools/{scale}/s2_plus.bigWig",
        ]

    def test_spikein_paths_use_the_spikein_directory(self, groupings):
        paths = get_condition_input_bigwigs(
            _Wildcards("treated"),
            pileup_method="bamnado",
            spikein_method="orlando",
            output_dir="out",
            sample_groupings=groupings,
            strand="minus",
        )
        assert paths == [
            "out/bigwigs/bamnado/spikein/orlando/s1_minus.bigWig",
            "out/bigwigs/bamnado/spikein/orlando/s2_minus.bigWig",
        ]

    def test_no_unexpanded_placeholders_remain(self, groupings):
        """Guards the specific failure mode: a non-f-string segment leaving
        '{sample}'/'{suffix}' in the template."""
        for kwargs in ({}, {"strand": "plus"}, {"spikein_method": "orlando"}):
            paths = get_condition_input_bigwigs(
                _Wildcards("treated"),
                pileup_method="deeptools",
                output_dir="out",
                sample_groupings=groupings,
                **kwargs,
            )
            assert paths, f"no paths returned for {kwargs}"
            for p in paths:
                assert "{" not in p and "}" not in p, f"unexpanded placeholder in {p}"
