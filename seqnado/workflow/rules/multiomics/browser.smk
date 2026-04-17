
# Define input and output files so that seqnado isnt required in the rule

plot_files=OUTPUT.track_plots


rule plotnado_multiomics_deeptools:
    input:
        rules.gather_bigwigs.output.bw_dir,
    output:
        plots=plot_files,
        template=OUTPUT_DIR + "/multiomics/track_plots/template.toml",
    params:
        assay=CONFIG.assay.value,
        genes=str(CONFIG.genome.genes) if CONFIG.assay_config.plot_with_plotnado and CONFIG.genome.genes else None,
        outdir=OUTPUT_DIR + "/multiomics/track_plots/",
        regions=CONFIG.assay_config.plotting.coordinates if CONFIG.assay_config.plotting else None,
        plotting_format=CONFIG.assay_config.plotting.file_format if CONFIG.assay_config.plotting else None,
    resources:
        mem="1.5GB",
         runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
    container: "library://asmith151/plotnado/plotnado:latest"
    log: OUTPUT_DIR + "/multiomics/logs/visualise/plotnado.log",
    benchmark: OUTPUT_DIR + ".benchmark/multiomics/visualise/plotnado.tsv",
    message: "Generating genome browser visualisations with Plotnado"
    script:
        "../../scripts/run_plotnado.py"
