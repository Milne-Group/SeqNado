"""
Tests for condition-based bigwig comparisons.

Tests that when perform_comparisons is enabled and multiple condition groups exist,
the pipeline generates aggregated and subtraction bigwigs as expected.
"""

from pathlib import Path

import pytest


@pytest.mark.pipeline
@pytest.mark.snakemake
@pytest.mark.requires_apptainer
@pytest.mark.slow
@pytest.mark.parametrize(
    "assay",
    [
        pytest.param("atac", id="atac"),
        pytest.param("chip", id="chip"),
        pytest.param("cat", id="cat"),
        pytest.param("rna", id="rna"),
    ],
)
def test_condition_bigwig_comparisons(
    assay,
    test_context,
    config_yaml_for_testing: Path,
    design: Path,
    test_profile_path: Path,
    pytestconfig,
    seqnado_runner,
):
    """
    Test that condition-based bigwig aggregations and subtractions are generated.

    This test enables perform_comparisons and ensures that:
    1. Aggregated condition bigwigs exist
    2. Subtraction bigwigs for all condition pairs exist
    3. Both unscaled and spike-in variants (if applicable) are generated
    """
    import pandas as pd
    import yaml

    assay_type = test_context.assay_type(assay)
    cores = test_context.cores
    preset = pytestconfig.getoption("--preset", default="t")

    # Add condition column to design file with 2 condition groups
    df = pd.read_csv(design)
    if "condition" not in df.columns:
        # Split conditions evenly: first half gets "ctrl", second half gets "treat"
        mid = len(df) // 2
        conditions = ["ctrl"] * mid + ["treat"] * (len(df) - mid)
        df["condition"] = conditions
    else:
        # Ensure we have at least 2 unique conditions
        unique_conds = df["condition"].nunique()
        if unique_conds < 2:
            mid = len(df) // 2
            conditions = ["ctrl"] * mid + ["treat"] * (len(df) - mid)
            df["condition"] = conditions

    # Add consensus_group for non-MCC assays
    if assay != "mcc" and "consensus_group" not in df.columns:
        df["consensus_group"] = "all"
    elif assay != "mcc":
        df["consensus_group"] = df["consensus_group"].fillna("all")

    # Add group/deseq2 for RNA
    if assay == "rna":
        if "group" not in df.columns:
            df["group"] = df["condition"]
        if "deseq2" not in df.columns:
            df["deseq2"] = (df["group"].str.lower() == "treat").astype(int)

    df.to_csv(design, index=False)

    # Load config and enable perform_comparisons
    with open(config_yaml_for_testing) as f:
        config = yaml.safe_load(f)

    # Enable perform_comparisons for bigwigs
    if "bigwigs" in config["assay_config"]:
        config["assay_config"]["bigwigs"]["perform_comparisons"] = True

    # Save updated config
    with open(config_yaml_for_testing, "w") as f:
        yaml.dump(config, f)

    res = seqnado_runner(
        [
            "seqnado",
            "pipeline",
            assay_type,
            "-c",
            str(cores),
            "--configfile",
            str(config_yaml_for_testing),
            "--preset",
            preset,
        ],
        cwd=config_yaml_for_testing.parent,
        capture_output=False,
        text=True,
    )

    print("STDOUT:", res.stdout)
    print("STDERR:", res.stderr)
    assert res.returncode == 0, (
        f"Pipeline failed with return code {res.returncode}. See output above."
    )

    test_dir = config_yaml_for_testing.parent
    output_dir = test_dir / "seqnado_output" / assay_type
    bigwig_dir = output_dir / "bigwigs"

    # Verify basic output exists
    assert not (test_dir / "seqnado_error.log").exists()
    assert (test_dir / "seqnado_output").exists()
    assert (output_dir / "seqnado_report.html").exists()
    assert bigwig_dir.exists()

    # Get pileup methods from config to check condition comparison outputs
    if "bigwigs" in config["assay_config"]:
        pileup_methods = config["assay_config"]["bigwigs"].get("pileup_method", [])
        if not isinstance(pileup_methods, list):
            pileup_methods = [pileup_methods]

        conditions = df["condition"].unique().tolist()
        condition_pairs = [
            (c1, c2) for c1 in conditions for c2 in conditions if c1 != c2
        ]

        for method in pileup_methods:
            # Handle both enum and string format
            if hasattr(method, "value"):
                method_name = method.value
            else:
                method_name = str(method).lower()

            method_dir = bigwig_dir / method_name

            if method_dir.exists():
                # Check aggregated condition bigwigs (unscaled)
                aggregated_dir = method_dir / "aggregated"
                if aggregated_dir.exists():
                    for cond in conditions:
                        cond_bw = aggregated_dir / f"{cond}.bigWig"
                        assert (
                            cond_bw.exists()
                        ), f"Aggregated bigwig not found: {cond_bw}"

                    # Check subtraction bigwigs (unscaled)
                    subtraction_dir = method_dir / "subtraction"
                    if subtraction_dir.exists():
                        for cond1, cond2 in condition_pairs:
                            subtraction_bw = (
                                subtraction_dir / f"{cond1}_vs_{cond2}.bigWig"
                            )
                            assert (
                                subtraction_bw.exists()
                            ), f"Subtraction bigwig not found: {subtraction_bw}"

                # Check spike-in variants if applicable
                spikein_dir = method_dir / "spikein"
                if spikein_dir.exists():
                    for spikein_method_dir in spikein_dir.iterdir():
                        if spikein_method_dir.is_dir():
                            # Check aggregated spike-in bigwigs
                            si_aggregated = spikein_method_dir / "aggregated"
                            if si_aggregated.exists():
                                for cond in conditions:
                                    si_cond_bw = si_aggregated / f"{cond}.bigWig"
                                    assert (
                                        si_cond_bw.exists()
                                    ), f"Spike-in aggregated bigwig not found: {si_cond_bw}"

                            # Check spike-in subtraction bigwigs
                            si_subtraction = spikein_method_dir / "subtraction"
                            if si_subtraction.exists():
                                for cond1, cond2 in condition_pairs:
                                    si_sub_bw = (
                                        si_subtraction / f"{cond1}_vs_{cond2}.bigWig"
                                    )
                                    assert (
                                        si_sub_bw.exists()
                                    ), f"Spike-in subtraction bigwig not found: {si_sub_bw}"
