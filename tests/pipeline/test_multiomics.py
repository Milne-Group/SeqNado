import subprocess
from pathlib import Path

import pytest


@pytest.mark.pipeline
@pytest.mark.snakemake
@pytest.mark.requires_apptainer
@pytest.mark.slow
@pytest.mark.parametrize(
    "multiomics",
    [
        [
            "atac",
            "chip",
            "rna",
        ]
    ],
)
def test_multiomics(
    multiomics: list[str],
    multiomics_configs: dict[str, dict[str, Path]],
    multiomics_run_directory: Path,
    cores: int,
):
    """Test running multiple assays in a single seqnado project using Snakefile_multi.

    This test uses the new SeqNado multiomics configuration design which includes:
    - Individual config_<assay>.yaml files for each assay (created via seqnado config)
    - A config_multiomics.yaml file specifying multiomics-specific settings
    - Automatic discovery of assay configs by Snakefile_multi

    The test verifies that:
    1. Individual assay configs are properly created by the fixture (using seqnado config)
    2. The multiomics config is created
    3. The multiomics pipeline runs successfully
    4. All expected outputs are generated (individual assay outputs + multiomics outputs)

    Args:
        multiomics: List of assay names to run together
        multiomics_configs: Dict mapping assay names to their config and metadata paths
                            (configs created using seqnado config command)
        multiomics_run_directory: Path to the run directory for the Multiomic test
        cores: Number of cores to use for the pipeline
    """
    # Verify that config_multiomics.yaml was created by the fixture
    multiomics_config = multiomics_run_directory / "config_multiomics.yaml"
    assert multiomics_config.exists(), (
        f"config_multiomics.yaml not found at {multiomics_config}"
    )

    # Verify that individual assay configs exist (created via seqnado config)
    for assay in multiomics:
        assay_config = multiomics_run_directory / f"config_{assay}.yaml"
        assert assay_config.exists(), (
            f"config_{assay}.yaml not found (should be created by seqnado config)"
        )

        assay_metadata = multiomics_run_directory / f"metadata_{assay}.csv"
        assert assay_metadata.exists(), (
            f"metadata_{assay}.csv not found (should be created by seqnado design)"
        )

    # Run the multiomics pipeline
    res = subprocess.run(
        [
            "seqnado",
            "pipeline",
            "-c",
            str(cores),
            "--preset",
            "t",
        ],
        cwd=multiomics_run_directory,
        capture_output=False,
        text=True,
    )

    assert res.returncode == 0, (
        f"Pipeline failed with return code {res.returncode}. See output above."
    )

    # Verify individual assay outputs
    for assay in multiomics:
        assay_output_dir = multiomics_run_directory / f"seqnado_output/{assay}"
        assert assay_output_dir.exists(), (
            f"Output directory not found for {assay}: {assay_output_dir}"
        )

        seqnado_report = assay_output_dir / "seqnado_report.html"
        assert seqnado_report.exists(), (
            f"No seqnado_report.html file found for {assay}: {seqnado_report}"
        )

    # Verify multiomics-specific outputs (based on MultiomicsOutput class)
    output_dir = multiomics_run_directory / "seqnado_output"

    # Summary report
    summary_report = output_dir / "multiomics_summary.txt"
    assert summary_report.exists(), (
        f"multiomics_summary.txt not found at {summary_report}"
    )

    # Heatmap outputs
    heatmap_pdf = output_dir / "multiomics" / "heatmap" / "heatmap.pdf"
    assert heatmap_pdf.exists(), (
        f"Multiomics heatmap not found at {heatmap_pdf}"
    )

    metaplot_pdf = output_dir / "multiomics" / "heatmap" / "metaplot.pdf"
    assert metaplot_pdf.exists(), (
        f"Multiomics metaplot not found at {metaplot_pdf}"
    )

    # ML dataset (using regions_bed mode with hg38_genes.bed)
    dataset_h5ad = output_dir / "multiomics" / "dataset" / "dataset_regions.h5ad"
    assert dataset_h5ad.exists(), (
        f"Multiomics dataset not found at {dataset_h5ad}"
    )
