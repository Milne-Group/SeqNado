from seqnado.workflow.helpers.common import define_time_requested
from seqnado import DataScalingTechnique, PileupMethod

# Genome browser plot generation from bigWig files using Plotnado.
# One rule per (pileup method) × (scale) × (merged) combination.
# Rules are only defined when the corresponding PlotFiles are registered in OUTPUT
# (i.e. when the user has configured plotnado for that combination).
# Compatible combinations:
#   deeptools  : individual → unscaled/csaw/spikein ; merged → unscaled/csaw/spikein
#   bamnado    : individual → unscaled/csaw/spikein ; merged → unscaled/csaw/spikein
#   homer      : individual → unscaled              ; merged → unscaled


# ─── deeptools · individual · unscaled (base rule) ───────────────────────────

rule generate_plotnado_visualisation:
    input:
        data=OUTPUT.select_bigwig_subtype(
            method=PileupMethod.DEEPTOOLS,
            scale=DataScalingTechnique.UNSCALED,
            ip_only=True,
        ) + OUTPUT.peak_files,
    output:
        plots=OUTPUT.select_genome_browser_plots(DataScalingTechnique.UNSCALED),
        template=OUTPUT_DIR + "/genome_browser_plots/deeptools/unscaled/template.toml",
    params:
        assay=str(CONFIG.assay) if hasattr(CONFIG, 'assay') else None,
        genes=str(CONFIG.genome.genes) if CONFIG.assay_config.plot_with_plotnado and hasattr(CONFIG.genome, 'genes') else None,
        regions=str(CONFIG.assay_config.plotting.coordinates) if CONFIG.assay_config.plotting and hasattr(CONFIG.assay_config.plotting, 'coordinates') else None,
        plotting_format=str(CONFIG.assay_config.plotting.file_format) if CONFIG.assay_config.plotting and hasattr(CONFIG.assay_config.plotting, 'file_format') else None,
    resources:
        mem="1.5GB",
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
    container: "library://asmith151/plotnado/plotnado:latest"
    log: OUTPUT_DIR + "/logs/visualise/plotnado_deeptools_unscaled.log",
    benchmark: OUTPUT_DIR + "/.benchmark/visualise/plotnado_deeptools_unscaled.tsv",
    message: "Generating deeptools unscaled genome browser visualisations with Plotnado"
    script:
        "../../scripts/run_plotnado.py"


# ─── deeptools · individual · csaw ───────────────────────────────────────────

if OUTPUT.select_genome_browser_plots(DataScalingTechnique.CSAW):
    use rule generate_plotnado_visualisation as generate_plotnado_visualisation_csaw with:
        input:
            data=OUTPUT.select_bigwig_subtype(
                method=PileupMethod.DEEPTOOLS,
                scale=DataScalingTechnique.CSAW,
                ip_only=True,
            ) + OUTPUT.peak_files,
        output:
            plots=OUTPUT.select_genome_browser_plots(DataScalingTechnique.CSAW),
            template=OUTPUT_DIR + "/genome_browser_plots/deeptools/csaw/template.toml",
        log: OUTPUT_DIR + "/logs/visualise/plotnado_deeptools_csaw.log",
        benchmark: OUTPUT_DIR + "/.benchmark/visualise/plotnado_deeptools_csaw.tsv",
        message: "Generating deeptools CSAW-scaled genome browser visualisations with Plotnado"


# ─── deeptools · individual · spikein ────────────────────────────────────────

if OUTPUT.select_genome_browser_plots(DataScalingTechnique.SPIKEIN):
    use rule generate_plotnado_visualisation as generate_plotnado_visualisation_spikein with:
        input:
            data=OUTPUT.select_bigwig_subtype(
                method=PileupMethod.DEEPTOOLS,
                scale=DataScalingTechnique.SPIKEIN,
                ip_only=True,
            ) + OUTPUT.peak_files,
        output:
            plots=OUTPUT.select_genome_browser_plots(DataScalingTechnique.SPIKEIN),
            template=OUTPUT_DIR + "/genome_browser_plots/deeptools/spikein/template.toml",
        log: OUTPUT_DIR + "/logs/visualise/plotnado_deeptools_spikein.log",
        benchmark: OUTPUT_DIR + "/.benchmark/visualise/plotnado_deeptools_spikein.tsv",
        message: "Generating deeptools spike-in normalized genome browser visualisations with Plotnado"


# ─── deeptools · merged · unscaled ───────────────────────────────────────────

if OUTPUT.select_genome_browser_plots(DataScalingTechnique.UNSCALED, is_merged=True):
    use rule generate_plotnado_visualisation as generate_plotnado_visualisation_merged with:
        input:
            data=OUTPUT.select_bigwig_subtype(
                method=PileupMethod.DEEPTOOLS,
                scale=DataScalingTechnique.UNSCALED,
                is_merged=True,
            ) + OUTPUT.peak_files,
        output:
            plots=OUTPUT.select_genome_browser_plots(DataScalingTechnique.UNSCALED, is_merged=True),
            template=OUTPUT_DIR + "/genome_browser_plots/merged/deeptools/unscaled/template.toml",
        log: OUTPUT_DIR + "/logs/visualise/plotnado_merged_deeptools_unscaled.log",
        benchmark: OUTPUT_DIR + "/.benchmark/visualise/plotnado_merged_deeptools_unscaled.tsv",
        message: "Generating deeptools merged unscaled genome browser visualisations with Plotnado"


# ─── deeptools · merged · csaw ───────────────────────────────────────────────

if OUTPUT.select_genome_browser_plots(DataScalingTechnique.CSAW, is_merged=True):
    use rule generate_plotnado_visualisation as generate_plotnado_visualisation_merged_csaw with:
        input:
            data=OUTPUT.select_bigwig_subtype(
                method=PileupMethod.DEEPTOOLS,
                scale=DataScalingTechnique.CSAW,
                is_merged=True,
            ) + OUTPUT.peak_files,
        output:
            plots=OUTPUT.select_genome_browser_plots(DataScalingTechnique.CSAW, is_merged=True),
            template=OUTPUT_DIR + "/genome_browser_plots/merged/deeptools/csaw/template.toml",
        log: OUTPUT_DIR + "/logs/visualise/plotnado_merged_deeptools_csaw.log",
        benchmark: OUTPUT_DIR + "/.benchmark/visualise/plotnado_merged_deeptools_csaw.tsv",
        message: "Generating deeptools merged CSAW-scaled genome browser visualisations with Plotnado"


# ─── deeptools · merged · spikein ────────────────────────────────────────────

if OUTPUT.select_genome_browser_plots(DataScalingTechnique.SPIKEIN, is_merged=True):
    use rule generate_plotnado_visualisation as generate_plotnado_visualisation_merged_spikein with:
        input:
            data=OUTPUT.select_bigwig_subtype(
                method=PileupMethod.DEEPTOOLS,
                scale=DataScalingTechnique.SPIKEIN,
                is_merged=True,
            ) + OUTPUT.peak_files,
        output:
            plots=OUTPUT.select_genome_browser_plots(DataScalingTechnique.SPIKEIN, is_merged=True),
            template=OUTPUT_DIR + "/genome_browser_plots/merged/deeptools/spikein/template.toml",
        log: OUTPUT_DIR + "/logs/visualise/plotnado_merged_deeptools_spikein.log",
        benchmark: OUTPUT_DIR + "/.benchmark/visualise/plotnado_merged_deeptools_spikein.tsv",
        message: "Generating deeptools merged spike-in normalized genome browser visualisations with Plotnado"


# ─── bamnado · individual · unscaled ─────────────────────────────────────────

if OUTPUT.select_genome_browser_plots(DataScalingTechnique.UNSCALED, method=PileupMethod.BAMNADO):
    use rule generate_plotnado_visualisation as generate_plotnado_visualisation_bamnado with:
        input:
            data=OUTPUT.select_bigwig_subtype(
                method=PileupMethod.BAMNADO,
                scale=DataScalingTechnique.UNSCALED,
                ip_only=True,
            ) + OUTPUT.peak_files,
        output:
            plots=OUTPUT.select_genome_browser_plots(DataScalingTechnique.UNSCALED, method=PileupMethod.BAMNADO),
            template=OUTPUT_DIR + "/genome_browser_plots/bamnado/unscaled/template.toml",
        log: OUTPUT_DIR + "/logs/visualise/plotnado_bamnado_unscaled.log",
        benchmark: OUTPUT_DIR + "/.benchmark/visualise/plotnado_bamnado_unscaled.tsv",
        message: "Generating bamnado unscaled genome browser visualisations with Plotnado"


# ─── bamnado · merged · unscaled ─────────────────────────────────────────────

if OUTPUT.select_genome_browser_plots(DataScalingTechnique.UNSCALED, method=PileupMethod.BAMNADO, is_merged=True):
    use rule generate_plotnado_visualisation as generate_plotnado_visualisation_bamnado_merged with:
        input:
            data=OUTPUT.select_bigwig_subtype(
                method=PileupMethod.BAMNADO,
                scale=DataScalingTechnique.UNSCALED,
                is_merged=True,
            ) + OUTPUT.peak_files,
        output:
            plots=OUTPUT.select_genome_browser_plots(DataScalingTechnique.UNSCALED, method=PileupMethod.BAMNADO, is_merged=True),
            template=OUTPUT_DIR + "/genome_browser_plots/merged/bamnado/unscaled/template.toml",
        log: OUTPUT_DIR + "/logs/visualise/plotnado_merged_bamnado_unscaled.log",
        benchmark: OUTPUT_DIR + "/.benchmark/visualise/plotnado_merged_bamnado_unscaled.tsv",
        message: "Generating bamnado merged unscaled genome browser visualisations with Plotnado"


# ─── bamnado · merged · csaw ─────────────────────────────────────────────────

if OUTPUT.select_genome_browser_plots(DataScalingTechnique.CSAW, method=PileupMethod.BAMNADO, is_merged=True):
    use rule generate_plotnado_visualisation as generate_plotnado_visualisation_bamnado_merged_csaw with:
        input:
            data=OUTPUT.select_bigwig_subtype(
                method=PileupMethod.BAMNADO,
                scale=DataScalingTechnique.CSAW,
                is_merged=True,
            ) + OUTPUT.peak_files,
        output:
            plots=OUTPUT.select_genome_browser_plots(DataScalingTechnique.CSAW, method=PileupMethod.BAMNADO, is_merged=True),
            template=OUTPUT_DIR + "/genome_browser_plots/merged/bamnado/csaw/template.toml",
        log: OUTPUT_DIR + "/logs/visualise/plotnado_merged_bamnado_csaw.log",
        benchmark: OUTPUT_DIR + "/.benchmark/visualise/plotnado_merged_bamnado_csaw.tsv",
        message: "Generating bamnado merged CSAW-scaled genome browser visualisations with Plotnado"


# ─── bamnado · merged · spikein ──────────────────────────────────────────────

if OUTPUT.select_genome_browser_plots(DataScalingTechnique.SPIKEIN, method=PileupMethod.BAMNADO, is_merged=True):
    use rule generate_plotnado_visualisation as generate_plotnado_visualisation_bamnado_merged_spikein with:
        input:
            data=OUTPUT.select_bigwig_subtype(
                method=PileupMethod.BAMNADO,
                scale=DataScalingTechnique.SPIKEIN,
                is_merged=True,
            ) + OUTPUT.peak_files,
        output:
            plots=OUTPUT.select_genome_browser_plots(DataScalingTechnique.SPIKEIN, method=PileupMethod.BAMNADO, is_merged=True),
            template=OUTPUT_DIR + "/genome_browser_plots/merged/bamnado/spikein/template.toml",
        log: OUTPUT_DIR + "/logs/visualise/plotnado_merged_bamnado_spikein.log",
        benchmark: OUTPUT_DIR + "/.benchmark/visualise/plotnado_merged_bamnado_spikein.tsv",
        message: "Generating bamnado merged spike-in normalized genome browser visualisations with Plotnado"


# ─── homer · individual · unscaled ───────────────────────────────────────────

if OUTPUT.select_genome_browser_plots(DataScalingTechnique.UNSCALED, method=PileupMethod.HOMER):
    use rule generate_plotnado_visualisation as generate_plotnado_visualisation_homer with:
        input:
            data=OUTPUT.select_bigwig_subtype(
                method=PileupMethod.HOMER,
                scale=DataScalingTechnique.UNSCALED,
                ip_only=True,
            ) + OUTPUT.peak_files,
        output:
            plots=OUTPUT.select_genome_browser_plots(DataScalingTechnique.UNSCALED, method=PileupMethod.HOMER),
            template=OUTPUT_DIR + "/genome_browser_plots/homer/unscaled/template.toml",
        log: OUTPUT_DIR + "/logs/visualise/plotnado_homer_unscaled.log",
        benchmark: OUTPUT_DIR + "/.benchmark/visualise/plotnado_homer_unscaled.tsv",
        message: "Generating homer unscaled genome browser visualisations with Plotnado"


# ─── homer · merged · unscaled ───────────────────────────────────────────────

if OUTPUT.select_genome_browser_plots(DataScalingTechnique.UNSCALED, method=PileupMethod.HOMER, is_merged=True):
    use rule generate_plotnado_visualisation as generate_plotnado_visualisation_homer_merged with:
        input:
            data=OUTPUT.select_bigwig_subtype(
                method=PileupMethod.HOMER,
                scale=DataScalingTechnique.UNSCALED,
                is_merged=True,
            ) + OUTPUT.peak_files,
        output:
            plots=OUTPUT.select_genome_browser_plots(DataScalingTechnique.UNSCALED, method=PileupMethod.HOMER, is_merged=True),
            template=OUTPUT_DIR + "/genome_browser_plots/merged/homer/unscaled/template.toml",
        log: OUTPUT_DIR + "/logs/visualise/plotnado_merged_homer_unscaled.log",
        benchmark: OUTPUT_DIR + "/.benchmark/visualise/plotnado_merged_homer_unscaled.tsv",
        message: "Generating homer merged unscaled genome browser visualisations with Plotnado"
