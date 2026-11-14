"""Unit tests for MultiomicsOutputBuilder."""

from pathlib import Path
import pytest
from seqnado import Assay, PileupMethod, DataScalingTechnique
from seqnado.outputs.core import MultiomicsOutputBuilder, SeqnadoOutputFiles

class TestMultiomicsOutputBuilder:
    """Tests for MultiomicsOutputBuilder."""

    def test_initialization(self):
        """Test initialization with default and custom values."""
        # Default
        builder = MultiomicsOutputBuilder()
        assert builder.output_dir == Path("seqnado_output")
        assert builder.assay_outputs == {}
        assert builder.file_collections == []

        # Custom
        custom_dir = Path("custom_output")
        assay_outputs = {Assay.ATAC: SeqnadoOutputFiles(files=["file1"])}
        builder = MultiomicsOutputBuilder(output_dir=custom_dir, assay_outputs=assay_outputs)
        assert builder.output_dir == custom_dir
        assert builder.assay_outputs == assay_outputs

    def test_add_assay_bigwigs(self):
        """Test adding assay bigwigs."""
        assay_outputs = {
            Assay.ATAC: SeqnadoOutputFiles(files=["atac.bigWig", "atac.txt"]),
            Assay.RNA: SeqnadoOutputFiles(files=["rna.bigWig", "rna.bam"]),
        }
        builder = MultiomicsOutputBuilder(assay_outputs=assay_outputs)
        
        builder.add_assay_bigwigs()
        
        assert len(builder.file_collections) == 1
        files = builder.file_collections[0].files
        assert len(files) == 2
        assert "atac.bigWig" in files
        assert "rna.bigWig" in files
        assert "atac.txt" not in files

    def test_add_assay_peaks(self):
        """Test adding assay peaks."""
        # Mocking peak files property which relies on suffix filtering
        assay_outputs = {
            Assay.ATAC: SeqnadoOutputFiles(files=["atac_peaks.bed", "atac.txt"]),
            Assay.CHIP: SeqnadoOutputFiles(files=["chip_peaks.bed", "chip.bam"]),
        }
        builder = MultiomicsOutputBuilder(assay_outputs=assay_outputs)
        
        builder.add_assay_peaks()
        
        assert len(builder.file_collections) == 1
        files = builder.file_collections[0].files
        # Note: SeqnadoOutputFiles.peak_files property filters for .bed
        assert "atac_peaks.bed" in files
        assert "chip_peaks.bed" in files

    def test_add_summary_report(self):
        """Test adding summary report."""
        builder = MultiomicsOutputBuilder(output_dir=Path("out"))
        builder.add_summary_report()
        
        assert len(builder.file_collections) == 1
        assert str(Path("out/multiomics_summary.txt")) in builder.file_collections[0].files

    def test_add_heatmap(self):
        """Test adding heatmap."""
        builder = MultiomicsOutputBuilder(output_dir=Path("out"))
        builder.add_heatmap()
        
        assert len(builder.file_collections) == 1
        assert str(Path("out/multiomics/heatmap/heatmap.pdf")) in builder.file_collections[0].files

    def test_add_metaplot(self):
        """Test adding metaplot."""
        builder = MultiomicsOutputBuilder(output_dir=Path("out"))
        builder.add_metaplot()
        
        assert len(builder.file_collections) == 1
        assert str(Path("out/multiomics/heatmap/metaplot.pdf")) in builder.file_collections[0].files

    def test_add_multiomics_dataset(self):
        """Test adding multiomics dataset."""
        builder = MultiomicsOutputBuilder(output_dir=Path("out"))
        builder.add_multiomics_dataset()
        
        assert len(builder.file_collections) == 1
        assert str(Path("out/dataset/dataset.zarr")) in builder.file_collections[0].files

    def test_add_assay_outputs(self):
        """Test adding all assay outputs."""
        assay_outputs = {
            Assay.ATAC: SeqnadoOutputFiles(files=["f1", "f2"]),
            Assay.RNA: SeqnadoOutputFiles(files=["f3"]),
        }
        builder = MultiomicsOutputBuilder(assay_outputs=assay_outputs)
        
        builder.add_assay_outputs()
        
        assert len(builder.file_collections) == 2
        all_files = []
        for col in builder.file_collections:
            all_files.extend(col.files)
        
        assert len(all_files) == 3
        assert "f1" in all_files
        assert "f2" in all_files
        assert "f3" in all_files

    def test_build(self):
        """Test building the final SeqnadoOutputFiles."""
        builder = MultiomicsOutputBuilder(output_dir=Path("out"))
        builder.add_summary_report()
        builder.add_heatmap()
        
        output = builder.build()
        
        assert isinstance(output, SeqnadoOutputFiles)
        assert len(output.files) == 2
        assert str(Path("out/multiomics_summary.txt")) in output.files
        assert str(Path("out/multiomics/heatmap/heatmap.pdf")) in output.files



def _plotting(coordinates, file_format="png"):
    from seqnado.config.configs import PlottingConfig

    return PlottingConfig(coordinates=str(coordinates), file_format=file_format)


@pytest.mark.unit
class TestMultiomicsBrowserPlots:
    """Tests for the combined multiomics genome-browser plot outputs."""

    @pytest.fixture
    def regions(self, tmp_path):
        bed = tmp_path / "regions.bed"
        bed.write_text("chr1\t100\t200\tlocusA\nchr2\t300\t400\tlocusB\n")
        return bed

    @pytest.fixture
    def assay_outputs(self):
        return {
            Assay.CHIP: SeqnadoOutputFiles(
                files=[
                    "out/chip/bigwigs/deeptools/scaled-per-sample/ip1.bigWig",
                    "out/chip/bigwigs/deeptools/scaled-per-sample/input1.bigWig",
                    "out/chip/bigwigs/bamnado/scaled-per-sample/ip1.bigWig",
                    "out/chip/bigwigs/deeptools/merged/scaled-per-sample/merged/grp.bigWig",
                ],
                ip_sample_names=["ip1"],
            ),
            # No deeptools bigWigs: falls back to the next pileup method.
            Assay.ATAC: SeqnadoOutputFiles(
                files=["out/atac/bigwigs/bamnado/scaled-per-sample/atac1.bigWig"],
            ),
        }

    def test_browser_tracks_pick_one_method_per_assay(self, assay_outputs):
        builder = MultiomicsOutputBuilder(assay_outputs=assay_outputs)

        assert builder.browser_tracks == [
            ("out/chip/bigwigs/deeptools/scaled-per-sample/ip1.bigWig", Assay.CHIP),
            ("out/atac/bigwigs/bamnado/scaled-per-sample/atac1.bigWig", Assay.ATAC),
        ]

    @pytest.mark.parametrize("as_dict", [False, True], ids=["model", "dict"])
    def test_add_browser_plots_one_plot_per_region(self, assay_outputs, regions, as_dict):
        """All assays share one figure per region (no per-assay plots)."""
        plotting = _plotting(regions, "svg")
        builder = MultiomicsOutputBuilder(output_dir=Path("out"), assay_outputs=assay_outputs)

        assert builder.add_browser_plots(plotting.model_dump() if as_dict else plotting) is True
        assert builder.file_collections[0].files == [
            str(Path("out/multiomics/track_plots/locusA.svg")),
            str(Path("out/multiomics/track_plots/locusB.svg")),
            str(Path("out/multiomics/track_plots/template.toml")),
        ]

    @pytest.mark.parametrize(
        "has_plotting,has_bigwigs",
        [(False, True), (True, False)],
        ids=["no-plotting-config", "no-bigwigs"],
    )
    def test_add_browser_plots_skipped(self, regions, has_plotting, has_bigwigs):
        files = ["out/chip/bigwigs/deeptools/scaled-per-sample/s.bigWig"] if has_bigwigs else []
        builder = MultiomicsOutputBuilder(
            assay_outputs={Assay.CHIP: SeqnadoOutputFiles(files=files)},
        )

        assert builder.add_browser_plots(_plotting(regions) if has_plotting else None) is False
        assert builder.file_collections == []


@pytest.mark.unit
def test_take_browser_plotting_moves_first_plotting_and_disables_per_assay():
    from seqnado.outputs.multiomics import take_browser_plotting

    raw = {
        Assay.ATAC: {"assay_config": {"plot_with_plotnado": False}},
        Assay.CHIP: {
            "assay_config": {
                "plot_with_plotnado": True,
                "plotting": {"coordinates": "chip.bed", "file_format": "png"},
            }
        },
        Assay.METH: {
            "assay_config": {"plotting": {"coordinates": "meth.bed"}},
        },
        Assay.SNP: {"assay_config": None},
    }

    assert take_browser_plotting(raw) == {"coordinates": "chip.bed", "file_format": "png"}
    for config in raw.values():
        assay_config = config["assay_config"] or {}
        assert "plotting" not in assay_config
        assert assay_config.get("plot_with_plotnado", False) is False


@pytest.mark.unit
def test_take_browser_plotting_none_when_no_assay_plots():
    from seqnado.outputs.multiomics import take_browser_plotting

    assert take_browser_plotting({Assay.ATAC: {"assay_config": {}}}) is None
