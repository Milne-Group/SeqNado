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

    if assay == 'mcc':
        import pandas as pd
        # Update the `condition` column in the design file for mcc test so we can test condition-based groupings
        df = pd.read_csv(design)
        df['condition'] = df['sample_id'].str.split('-').str[1].str.split('_').str[0]
        df.to_csv(design, index=False)
    elif assay in ('atac', 'chip-rx', 'rna', 'rna-rx'):
        import pandas as pd
        # For ATAC-seq, ChIP-rx, and RNA assays, we need to add a 'consensus_group' column to the design file
        # for merged peak calling (ATAC/ChIP-rx) and general metadata validation (RNA)
        df = pd.read_csv(design)
        if 'consensus_group' not in df.columns:
            # For RNA-rx, use condition/group for consensus grouping so replicates are merged by condition
            if assay in ('rna-rx',):
                if 'condition' in df.columns:
                    df['consensus_group'] = df['condition']
                elif 'group' in df.columns:
                    df['consensus_group'] = df['group']
                else:
                    df['consensus_group'] = df['sample_id'].str.split('-').str[2]
            else:
                df['consensus_group'] = 'all'
        else:
            # Fill any NaN values in consensus_group with 'all'
            df['consensus_group'] = df['consensus_group'].fillna('all')
        # For RNA-RX, add group and deseq2 columns if not present
        if assay in ('rna', 'rna-rx'):
            if 'group' not in df.columns:
                df['group'] = df['sample_id'].str.split('-').str[2]
            if 'deseq2' not in df.columns:
                df['deseq2'] = (df['group'].str.lower() == 'treated').astype(int)
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


