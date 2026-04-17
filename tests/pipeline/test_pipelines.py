from pathlib import Path

import pytest


@pytest.mark.pipeline
@pytest.mark.snakemake
@pytest.mark.requires_apptainer
@pytest.mark.slow
def test_pipeline(
    assay,
    test_context,
    config_yaml_for_testing: Path,
    design: Path,
    test_profile_path: Path,
    pytestconfig,
    seqnado_runner,
):
    assay_type = test_context.assay_type(assay)
    cores = test_context.cores
    preset = pytestconfig.getoption("--preset", default="t")

    import pandas as pd
    df = pd.read_csv(design)

    # Add condition and metadata columns based on assay type
    if assay == 'mcc':
        # Extract condition from sample_id
        df['condition'] = df['sample_id'].str.split('-').str[1].str.split('_').str[0]
    elif assay in ('atac', 'chip', 'chip-rx', 'rna', 'rna-rx'):
        # Add conditions by splitting samples evenly: first half ctrl, second half treat
        if 'condition' not in df.columns:
            mid = len(df) // 2
            conditions = ['ctrl'] * mid + ['treat'] * (len(df) - mid)
            df['condition'] = conditions

        # Add consensus_group for peak calling assays
        if assay in ('atac', 'chip-rx'):
            if 'consensus_group' not in df.columns:
                df['consensus_group'] = 'all'
            else:
                df['consensus_group'] = df['consensus_group'].fillna('all')

        # For RNA-rx, use condition for consensus grouping to merge replicates by condition
        if assay == 'rna-rx':
            if 'consensus_group' not in df.columns:
                df['consensus_group'] = df['condition']
            else:
                df['consensus_group'] = df['consensus_group'].fillna('all')

        # Add group and deseq2 columns for RNA assays
        if assay in ('rna', 'rna-rx'):
            if 'group' not in df.columns:
                df['group'] = df['condition']
            if 'deseq2' not in df.columns:
                df['deseq2'] = (df['group'].str.lower() == 'treat').astype(int)

    df.to_csv(design, index=False)

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
            preset
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
    assert not (test_dir / "seqnado_error.log").exists()
    assert (test_dir / "seqnado_output").exists()
    assert (test_dir / f"seqnado_output/{assay_type}/seqnado_report.html").exists()


