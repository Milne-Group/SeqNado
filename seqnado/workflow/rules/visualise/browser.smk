from seqnado.workflow.helpers.common import define_time_requested, define_memory_requested
from seqnado import DataScalingTechnique, PileupMethod, SpikeInMethod

_spikein_cfg = getattr(CONFIG.assay_config, "spikein", None)
_spikein_methods = _spikein_cfg.method if _spikein_cfg else []
_has_spikein_orlando = SpikeInMethod.ORLANDO in _spikein_methods
_has_spikein_withinput = SpikeInMethod.WITH_INPUT in _spikein_methods

# Genome browser plot generation from bigWig files using Plotnado.
# One rule per (pileup method) × (scale) × (merged) combination.
# Rules are only defined when the corresponding PlotFiles are registered in OUTPUT
# (i.e. when the user has configured plotnado for that combination).
# Compatible combinations:
#   deeptools  : individual → unscaled/csaw/spikein_orlando/spikein_withinput ; merged → unscaled/csaw/spikein_orlando/spikein_withinput
#   bamnado    : individual → unscaled/csaw/spikein_orlando/spikein_withinput ; merged → unscaled/csaw/spikein_orlando/spikein_withinput
#   homer      : individual → unscaled              ; merged → unscaled

rule index_individual_peaks:
    input:
        peaks=OUTPUT.peak_files,
    output:
        peaks_indexed=temp(expand("{peak}.gz", peak=OUTPUT.peak_files)),
        peaks_tbi=temp(expand("{peak}.gz.tbi", peak=OUTPUT.peak_files)),
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
    log: OUTPUT_DIR + "/logs/visualise/index_individual_peaks.log",
    benchmark: OUTPUT_DIR + "/.benchmark/visualise/index_individual_peaks.tsv",
    message: "Indexing individual peak files"
    shell:
        """
        for bed in {input.peaks}; do
            sort -k1,1 -k2,2n "$bed" | bgzip > "$bed.gz"
            tabix -p bed "$bed.gz"
        done
        """


rule index_genes_bed:
    input:
        genes=str(CONFIG.genome.genes) if hasattr(CONFIG.genome, 'genes') else [],
    output:
        bed_gz=temp(OUTPUT_DIR + "/resources/genes_indexed.bed.gz"),
        tbi=temp(OUTPUT_DIR + "/resources/genes_indexed.bed.gz.tbi"),
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
    log: OUTPUT_DIR + "/logs/visualise/index_genes_bed.log",
    benchmark: OUTPUT_DIR + "/.benchmark/visualise/index_genes_bed.tsv",
    message: "Indexing genes file for plotting"
    shell:
        """
        if [ -s {input.genes} ]; then
            sort -k1,1 -k2,2n {input.genes} | bgzip > {output.bed_gz}
            tabix -p bed {output.bed_gz}
        else
            touch {output.bed_gz}
            touch {output.tbi}
        fi
        """


# ─── deeptools · individual · unscaled (base rule) ───────────────────────────
rule plotnado_deeptools:
    input:
        data=OUTPUT.select_bigwig_subtype(
            method=PileupMethod.DEEPTOOLS,
            scale=DataScalingTechnique.UNSCALED,
            ip_only=True,
        ),
        peaks_indexed=expand("{peak}.gz", peak=OUTPUT.peak_files),
        peaks_tbi=expand("{peak}.gz.tbi", peak=OUTPUT.peak_files),
        genes_indexed=rules.index_genes_bed.output.bed_gz,
        genes_indexed_tbi=rules.index_genes_bed.output.tbi,
    output:
        plots=OUTPUT.select_genome_browser_plots(DataScalingTechnique.UNSCALED),
        template=OUTPUT_DIR + "/genome_browser_plots/deeptools/unscaled/template.toml",
    params:
        assay=str(CONFIG.assay) if hasattr(CONFIG, 'assay') else None,
        peak_files=expand("{peak}.gz", peak=OUTPUT.peak_files),
        genes=OUTPUT_DIR + "/resources/genes_indexed.bed.gz" if CONFIG.assay_config.plot_with_plotnado and hasattr(CONFIG.genome, 'genes') else None,
        regions=str(CONFIG.assay_config.plotting.coordinates) if CONFIG.assay_config.plotting and hasattr(CONFIG.assay_config.plotting, 'coordinates') else None,
        plotting_format=str(CONFIG.assay_config.plotting.file_format) if CONFIG.assay_config.plotting and hasattr(CONFIG.assay_config.plotting, 'file_format') else None,
        outdir=lambda wildcards, output: str(os.path.dirname(output.template)),
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
use rule plotnado_deeptools as plotnado_deeptools_csaw with:
    input:
        data=OUTPUT.select_bigwig_subtype(
            method=PileupMethod.DEEPTOOLS,
            scale=DataScalingTechnique.CSAW,
            ip_only=True,
        ),
        peaks_indexed=expand("{peak}.gz", peak=OUTPUT.peak_files),
        peaks_tbi=expand("{peak}.gz.tbi", peak=OUTPUT.peak_files),
        genes_indexed=rules.index_genes_bed.output.bed_gz,
        genes_indexed_tbi=rules.index_genes_bed.output.tbi,
    output:
        plots=OUTPUT.select_genome_browser_plots(DataScalingTechnique.CSAW),
        template=OUTPUT_DIR + "/genome_browser_plots/deeptools/csaw/template.toml",
    log: OUTPUT_DIR + "/logs/visualise/plotnado_deeptools_csaw.log",
    benchmark: OUTPUT_DIR + "/.benchmark/visualise/plotnado_deeptools_csaw.tsv",
    message: "Generating deeptools CSAW-scaled genome browser visualisations with Plotnado"


# ─── deeptools · individual · spikein · orlando ───────────────────────────────
use rule plotnado_deeptools as plotnado_deeptools_spikein_orlando with:
    input:
        data=OUTPUT.select_bigwig_subtype(
            method=PileupMethod.DEEPTOOLS,
            scale=DataScalingTechnique.SPIKEIN,
            spikein_method=SpikeInMethod.ORLANDO.value,
            ip_only=True,
        ),
        peaks_indexed=expand("{peak}.gz", peak=OUTPUT.peak_files),
        peaks_tbi=expand("{peak}.gz.tbi", peak=OUTPUT.peak_files),
        genes_indexed=rules.index_genes_bed.output.bed_gz,
        genes_indexed_tbi=rules.index_genes_bed.output.tbi,
    output:
        plots=OUTPUT.select_genome_browser_plots(DataScalingTechnique.SPIKEIN, spikein_method=SpikeInMethod.ORLANDO.value),
        template=OUTPUT_DIR + "/genome_browser_plots/deeptools/spikein/orlando/template.toml",
    log: OUTPUT_DIR + "/logs/visualise/plotnado_deeptools_spikein_orlando.log",
    benchmark: OUTPUT_DIR + "/.benchmark/visualise/plotnado_deeptools_spikein_orlando.tsv",
    message: "Generating deeptools orlando spike-in genome browser visualisations with Plotnado"


# ─── deeptools · individual · spikein · with_input ───────────────────────────
use rule plotnado_deeptools as plotnado_deeptools_spikein_withinput with:
    input:
        data=OUTPUT.select_bigwig_subtype(
            method=PileupMethod.DEEPTOOLS,
            scale=DataScalingTechnique.SPIKEIN,
            spikein_method=SpikeInMethod.WITH_INPUT.value,
            ip_only=True,
        ),
        peaks_indexed=expand("{peak}.gz", peak=OUTPUT.peak_files),
        peaks_tbi=expand("{peak}.gz.tbi", peak=OUTPUT.peak_files),
        genes_indexed=rules.index_genes_bed.output.bed_gz,
        genes_indexed_tbi=rules.index_genes_bed.output.tbi,
    output:
        plots=OUTPUT.select_genome_browser_plots(DataScalingTechnique.SPIKEIN, spikein_method=SpikeInMethod.WITH_INPUT.value),
        template=OUTPUT_DIR + "/genome_browser_plots/deeptools/spikein/with_input/template.toml",
    log: OUTPUT_DIR + "/logs/visualise/plotnado_deeptools_spikein_withinput.log",
    benchmark: OUTPUT_DIR + "/.benchmark/visualise/plotnado_deeptools_spikein_withinput.tsv",
    message: "Generating deeptools with_input spike-in genome browser visualisations with Plotnado"


# ─── deeptools · merged · unscaled ───────────────────────────────────────────
use rule plotnado_deeptools as plotnado_deeptools_merged with:
    input:
        data=OUTPUT.select_bigwig_subtype(
            method=PileupMethod.DEEPTOOLS,
            scale=DataScalingTechnique.UNSCALED,
            is_merged=True,
        ),
        peaks_indexed=expand("{peak}.gz", peak=OUTPUT.peak_files),
        peaks_tbi=expand("{peak}.gz.tbi", peak=OUTPUT.peak_files),
        genes_indexed=rules.index_genes_bed.output.bed_gz,
        genes_indexed_tbi=rules.index_genes_bed.output.tbi,
    output:
        plots=OUTPUT.select_genome_browser_plots(DataScalingTechnique.UNSCALED, is_merged=True),
        template=OUTPUT_DIR + "/genome_browser_plots/merged/deeptools/unscaled/template.toml",
    log: OUTPUT_DIR + "/logs/visualise/plotnado_merged_deeptools_unscaled.log",
    benchmark: OUTPUT_DIR + "/.benchmark/visualise/plotnado_merged_deeptools_unscaled.tsv",
    message: "Generating deeptools merged unscaled genome browser visualisations with Plotnado"


# ─── deeptools · merged · csaw ───────────────────────────────────────────────
use rule plotnado_deeptools as plotnado_deeptools_merged_csaw with:
    input:
        data=OUTPUT.select_bigwig_subtype(
            method=PileupMethod.DEEPTOOLS,
            scale=DataScalingTechnique.CSAW,
            is_merged=True,
        ),
        peaks_indexed=expand("{peak}.gz", peak=OUTPUT.peak_files),
        peaks_tbi=expand("{peak}.gz.tbi", peak=OUTPUT.peak_files),
        genes_indexed=rules.index_genes_bed.output.bed_gz,
        genes_indexed_tbi=rules.index_genes_bed.output.tbi,
    output:
        plots=OUTPUT.select_genome_browser_plots(DataScalingTechnique.CSAW, is_merged=True),
        template=OUTPUT_DIR + "/genome_browser_plots/merged/deeptools/csaw/template.toml",
    log: OUTPUT_DIR + "/logs/visualise/plotnado_merged_deeptools_csaw.log",
    benchmark: OUTPUT_DIR + "/.benchmark/visualise/plotnado_merged_deeptools_csaw.tsv",
    message: "Generating deeptools merged CSAW-scaled genome browser visualisations with Plotnado"


# ─── deeptools · merged · spikein · orlando ──────────────────────────────────
use rule plotnado_deeptools as plotnado_deeptools_merged_spikein_orlando with:
    input:
        data=OUTPUT.select_bigwig_subtype(
            method=PileupMethod.DEEPTOOLS,
            scale=DataScalingTechnique.SPIKEIN,
            spikein_method=SpikeInMethod.ORLANDO.value,
            is_merged=True,
        ),
        peaks_indexed=expand("{peak}.gz", peak=OUTPUT.peak_files),
        peaks_tbi=expand("{peak}.gz.tbi", peak=OUTPUT.peak_files),
        genes_indexed=rules.index_genes_bed.output.bed_gz,
        genes_indexed_tbi=rules.index_genes_bed.output.tbi,
    output:
        plots=OUTPUT.select_genome_browser_plots(DataScalingTechnique.SPIKEIN, is_merged=True, spikein_method=SpikeInMethod.ORLANDO.value),
        template=OUTPUT_DIR + "/genome_browser_plots/merged/deeptools/spikein/orlando/template.toml",
    log: OUTPUT_DIR + "/logs/visualise/plotnado_merged_deeptools_spikein_orlando.log",
    benchmark: OUTPUT_DIR + "/.benchmark/visualise/plotnado_merged_deeptools_spikein_orlando.tsv",
    message: "Generating deeptools merged orlando spike-in genome browser visualisations with Plotnado"


# ─── deeptools · merged · spikein · with_input ───────────────────────────────
use rule plotnado_deeptools as plotnado_deeptools_merged_spikein_withinput with:
    input:
        data=OUTPUT.select_bigwig_subtype(
            method=PileupMethod.DEEPTOOLS,
            scale=DataScalingTechnique.SPIKEIN,
            spikein_method=SpikeInMethod.WITH_INPUT.value,
            is_merged=True,
        ),
        peaks_indexed=expand("{peak}.gz", peak=OUTPUT.peak_files),
        peaks_tbi=expand("{peak}.gz.tbi", peak=OUTPUT.peak_files),
        genes_indexed=rules.index_genes_bed.output.bed_gz,
        genes_indexed_tbi=rules.index_genes_bed.output.tbi,
    output:
        plots=OUTPUT.select_genome_browser_plots(DataScalingTechnique.SPIKEIN, is_merged=True, spikein_method=SpikeInMethod.WITH_INPUT.value),
        template=OUTPUT_DIR + "/genome_browser_plots/merged/deeptools/spikein/with_input/template.toml",
    log: OUTPUT_DIR + "/logs/visualise/plotnado_merged_deeptools_spikein_withinput.log",
    benchmark: OUTPUT_DIR + "/.benchmark/visualise/plotnado_merged_deeptools_spikein_withinput.tsv",
    message: "Generating deeptools merged with_input spike-in genome browser visualisations with Plotnado"


# ─── bamnado · individual · unscaled ─────────────────────────────────────────
use rule plotnado_deeptools as plotnado_bamnado with:
    input:
        data=OUTPUT.select_bigwig_subtype(
            method=PileupMethod.BAMNADO,
            scale=DataScalingTechnique.UNSCALED,
            ip_only=True,
        ),
        peaks_indexed=expand("{peak}.gz", peak=OUTPUT.peak_files),
        peaks_tbi=expand("{peak}.gz.tbi", peak=OUTPUT.peak_files),
        genes_indexed=rules.index_genes_bed.output.bed_gz,
        genes_indexed_tbi=rules.index_genes_bed.output.tbi,
    output:
        plots=OUTPUT.select_genome_browser_plots(DataScalingTechnique.UNSCALED, method=PileupMethod.BAMNADO),
        template=OUTPUT_DIR + "/genome_browser_plots/bamnado/unscaled/template.toml",
    log: OUTPUT_DIR + "/logs/visualise/plotnado_bamnado_unscaled.log",
    benchmark: OUTPUT_DIR + "/.benchmark/visualise/plotnado_bamnado_unscaled.tsv",
    message: "Generating bamnado unscaled genome browser visualisations with Plotnado"


# ─── bamnado · merged · unscaled ─────────────────────────────────────────────
use rule plotnado_deeptools as plotnado_bamnado_merged with:
    input:
        data=OUTPUT.select_bigwig_subtype(
            method=PileupMethod.BAMNADO,
            scale=DataScalingTechnique.UNSCALED,
            is_merged=True,
        ),
        peaks_indexed=expand("{peak}.gz", peak=OUTPUT.peak_files),
        peaks_tbi=expand("{peak}.gz.tbi", peak=OUTPUT.peak_files),
        genes_indexed=rules.index_genes_bed.output.bed_gz,
        genes_indexed_tbi=rules.index_genes_bed.output.tbi,
    output:
        plots=OUTPUT.select_genome_browser_plots(DataScalingTechnique.UNSCALED, method=PileupMethod.BAMNADO, is_merged=True),
        template=OUTPUT_DIR + "/genome_browser_plots/merged/bamnado/unscaled/template.toml",
    log: OUTPUT_DIR + "/logs/visualise/plotnado_merged_bamnado_unscaled.log",
    benchmark: OUTPUT_DIR + "/.benchmark/visualise/plotnado_merged_bamnado_unscaled.tsv",
    message: "Generating bamnado merged unscaled genome browser visualisations with Plotnado"


# ─── bamnado · merged · csaw ─────────────────────────────────────────────────
use rule plotnado_deeptools as plotnado_bamnado_merged_csaw with:
    input:
        data=OUTPUT.select_bigwig_subtype(
            method=PileupMethod.BAMNADO,
            scale=DataScalingTechnique.CSAW,
            is_merged=True,
        ),
        peaks_indexed=expand("{peak}.gz", peak=OUTPUT.peak_files),
        peaks_tbi=expand("{peak}.gz.tbi", peak=OUTPUT.peak_files),
        genes_indexed=rules.index_genes_bed.output.bed_gz,
        genes_indexed_tbi=rules.index_genes_bed.output.tbi,
    output:
        plots=OUTPUT.select_genome_browser_plots(DataScalingTechnique.CSAW, method=PileupMethod.BAMNADO, is_merged=True),
        template=OUTPUT_DIR + "/genome_browser_plots/merged/bamnado/csaw/template.toml",
    log: OUTPUT_DIR + "/logs/visualise/plotnado_merged_bamnado_csaw.log",
    benchmark: OUTPUT_DIR + "/.benchmark/visualise/plotnado_merged_bamnado_csaw.tsv",
    message: "Generating bamnado merged CSAW-scaled genome browser visualisations with Plotnado"


# ─── bamnado · merged · spikein · orlando ────────────────────────────────────
use rule plotnado_deeptools as plotnado_bamnado_merged_spikein_orlando with:
    input:
        data=OUTPUT.select_bigwig_subtype(
            method=PileupMethod.BAMNADO,
            scale=DataScalingTechnique.SPIKEIN,
            spikein_method=SpikeInMethod.ORLANDO.value,
            is_merged=True,
        ),
        peaks_indexed=expand("{peak}.gz", peak=OUTPUT.peak_files),
        peaks_tbi=expand("{peak}.gz.tbi", peak=OUTPUT.peak_files),
        genes_indexed=rules.index_genes_bed.output.bed_gz,
        genes_indexed_tbi=rules.index_genes_bed.output.tbi,
    output:
        plots=OUTPUT.select_genome_browser_plots(DataScalingTechnique.SPIKEIN, method=PileupMethod.BAMNADO, is_merged=True, spikein_method=SpikeInMethod.ORLANDO.value),
        template=OUTPUT_DIR + "/genome_browser_plots/merged/bamnado/spikein/orlando/template.toml",
    log: OUTPUT_DIR + "/logs/visualise/plotnado_merged_bamnado_spikein_orlando.log",
    benchmark: OUTPUT_DIR + "/.benchmark/visualise/plotnado_merged_bamnado_spikein_orlando.tsv",
    message: "Generating bamnado merged orlando spike-in genome browser visualisations with Plotnado"


# ─── bamnado · merged · spikein · with_input ─────────────────────────────────
use rule plotnado_deeptools as plotnado_bamnado_merged_spikein_withinput with:
    input:
        data=OUTPUT.select_bigwig_subtype(
            method=PileupMethod.BAMNADO,
            scale=DataScalingTechnique.SPIKEIN,
            spikein_method=SpikeInMethod.WITH_INPUT.value,
            is_merged=True,
        ),
        peaks_indexed=expand("{peak}.gz", peak=OUTPUT.peak_files),
        peaks_tbi=expand("{peak}.gz.tbi", peak=OUTPUT.peak_files),
        genes_indexed=rules.index_genes_bed.output.bed_gz,
        genes_indexed_tbi=rules.index_genes_bed.output.tbi,
    output:
        plots=OUTPUT.select_genome_browser_plots(DataScalingTechnique.SPIKEIN, method=PileupMethod.BAMNADO, is_merged=True, spikein_method=SpikeInMethod.WITH_INPUT.value),
        template=OUTPUT_DIR + "/genome_browser_plots/merged/bamnado/spikein/with_input/template.toml",
    log: OUTPUT_DIR + "/logs/visualise/plotnado_merged_bamnado_spikein_withinput.log",
    benchmark: OUTPUT_DIR + "/.benchmark/visualise/plotnado_merged_bamnado_spikein_withinput.tsv",
    message: "Generating bamnado merged with_input spike-in genome browser visualisations with Plotnado"


# ─── homer · individual · unscaled ───────────────────────────────────────────
use rule plotnado_deeptools as plotnado_homer with:
    input:
        data=OUTPUT.select_bigwig_subtype(
            method=PileupMethod.HOMER,
            scale=DataScalingTechnique.UNSCALED,
            ip_only=True,
        ),
        peaks_indexed=expand("{peak}.gz", peak=OUTPUT.peak_files),
        peaks_tbi=expand("{peak}.gz.tbi", peak=OUTPUT.peak_files),
        genes_indexed=rules.index_genes_bed.output.bed_gz,
        genes_indexed_tbi=rules.index_genes_bed.output.tbi,
    output:
        plots=OUTPUT.select_genome_browser_plots(DataScalingTechnique.UNSCALED, method=PileupMethod.HOMER),
        template=OUTPUT_DIR + "/genome_browser_plots/homer/unscaled/template.toml",
    log: OUTPUT_DIR + "/logs/visualise/plotnado_homer_unscaled.log",
    benchmark: OUTPUT_DIR + "/.benchmark/visualise/plotnado_homer_unscaled.tsv",
    message: "Generating homer unscaled genome browser visualisations with Plotnado"


# ─── homer · merged · unscaled ───────────────────────────────────────────────
use rule plotnado_deeptools as plotnado_homer_merged with:
    input:
        data=OUTPUT.select_bigwig_subtype(
            method=PileupMethod.HOMER,
            scale=DataScalingTechnique.UNSCALED,
            is_merged=True,
        ),
        peaks_indexed=expand("{peak}.gz", peak=OUTPUT.peak_files),
        peaks_tbi=expand("{peak}.gz.tbi", peak=OUTPUT.peak_files),
        genes_indexed=rules.index_genes_bed.output.bed_gz,
        genes_indexed_tbi=rules.index_genes_bed.output.tbi,
    output:
        plots=OUTPUT.select_genome_browser_plots(DataScalingTechnique.UNSCALED, method=PileupMethod.HOMER, is_merged=True),
        template=OUTPUT_DIR + "/genome_browser_plots/merged/homer/unscaled/template.toml",
    log: OUTPUT_DIR + "/logs/visualise/plotnado_merged_homer_unscaled.log",
    benchmark: OUTPUT_DIR + "/.benchmark/visualise/plotnado_merged_homer_unscaled.tsv",
    message: "Generating homer merged unscaled genome browser visualisations with Plotnado"
