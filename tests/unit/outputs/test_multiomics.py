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
        assert str(Path("out/multiomics/dataset.zarr")) in builder.file_collections[0].files

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
