from seqnado import Assay
from seqnado.workflow.helpers.common import define_time_requested

# Combined genome-browser plots: one figure per region holding every assay's
# tracks (one per-sample bigWig per IP sample per assay, see
# MultiomicsOutputBuilder.browser_tracks). Regions and file format come from
# BROWSER_PLOTTING (Snakefile_multi); per-assay plots are switched off.
BROWSER_TRACKS = multiomics_builder.browser_tracks

# Genes: reuse the sorted, tabix-indexed genes BED that an assay module builds
# (index_genes_bed lives in rules/visualise/browser.smk, which only these
# assays include).
_ASSAYS_WITH_GENES_INDEX = (Assay.ATAC, Assay.CAT, Assay.CHIP, Assay.RNA, Assay.METH)
BROWSER_GENES = next(
    (
        f"{OUTPUT_DIR}/{assay.clean_name}/resources/genes_indexed.bed.gz"
        for assay in ASSAYS
        if assay in _ASSAYS_WITH_GENES_INDEX and CONFIGS_PER_ASSAY[assay].genome.genes
    ),
    None,
)


rule plotnado_multiomics:
    input:
        data=[bigwig for bigwig, _ in BROWSER_TRACKS],
        genes=[BROWSER_GENES] if BROWSER_GENES else [],
    output:
        plots=multiomics_builder.browser_plot_files(BROWSER_PLOTTING),
        template=str(multiomics_builder.browser_plot_dir / "template.toml"),
    params:
        assay="MULTIOMICS",
        track_assays=[assay.value for _, assay in BROWSER_TRACKS],
        genes=BROWSER_GENES,
        outdir=str(multiomics_builder.browser_plot_dir),
        regions=BROWSER_PLOTTING.coordinates,
        plotting_format=BROWSER_PLOTTING.file_format,
    resources:
        mem="1.5GB",
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=1),
    container: CONTAINERS.PLOTNADO
    log: OUTPUT_DIR + "/multiomics/logs/visualise/plotnado.log",
    benchmark: OUTPUT_DIR + "/.benchmark/multiomics/visualise/plotnado.tsv",
    message: "Generating combined multiomics genome browser plots with PlotNado"
    script:
        "../../scripts/run_plotnado.py"
