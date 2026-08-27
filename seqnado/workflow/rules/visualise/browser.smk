from seqnado.workflow.helpers.common import define_time_requested, define_memory_requested
from seqnado import Assay, DataScalingTechnique, PileupMethod, GroupScalingMethod, SpikeInMethod

# Browser tracks for the per-group technique are generated for one bamnado scaling
# method only — additional configured group_scaling_methods still produce bigwigs,
# just not dedicated browser tracks. OUTPUT picks which one, so the rules and the
# expected-file list cannot drift apart.
_DEFAULT_SCALING_METHOD = OUTPUT.plot_scaling_method.value

_spikein_cfg = CONFIG.assay_config.spikein
_spikein_methods = _spikein_cfg.method if _spikein_cfg else []
_has_spikein_orlando = SpikeInMethod.ORLANDO in _spikein_methods
_has_spikein_withinput = SpikeInMethod.WITH_INPUT in _spikein_methods

# Genome browser plot generation from bigWig files using Plotnado.
# One rule per (pileup method) × (scale) × (merged) combination.
# Rules are only defined when the corresponding PlotFiles are registered in OUTPUT
# (i.e. when the user has configured plotnado for that combination).
# Compatible combinations:
#   deeptools    : individual → per_sample/per_group/spikein_orlando/spikein_withinput ; merged → per_sample/per_group/spikein_orlando/spikein_withinput
#   bamnado      : individual → per_sample/spikein_orlando/spikein_withinput      ; merged → per_sample/per_group/spikein_orlando/spikein_withinput
#   (per_group is genomics-only — excluded entirely for RNA-seq)
#   homer        : individual → per_sample              ; merged → per_sample
#   methyldackel : individual → per_sample (METH assay only, reads from bigwigs/taps or bigwigs/wgbs)

rule index_individual_peaks:
    input:
        peaks=OUTPUT.peak_files,
    output:
        peaks_indexed=expand("{peak}.gz", peak=OUTPUT.peak_files),
        peaks_tbi=expand("{peak}.gz.tbi", peak=OUTPUT.peak_files),
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
        genes=str(CONFIG.genome.genes) if CONFIG.genome.genes else [],
    output:
        bed_gz=OUTPUT_DIR + "/resources/genes_indexed.bed.gz",
        tbi=OUTPUT_DIR + "/resources/genes_indexed.bed.gz.tbi",
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


# ─── deeptools · individual · per_sample (base rule) ───────────────────────────
rule plotnado_deeptools:
    input:
        data=OUTPUT.select_bigwig_subtype(
            method=PileupMethod.DEEPTOOLS,
            scale=DataScalingTechnique.PER_SAMPLE,
            ip_only=True,
        ),
        peaks_indexed=expand("{peak}.gz", peak=OUTPUT.peak_files),
        peaks_tbi=expand("{peak}.gz.tbi", peak=OUTPUT.peak_files),
        genes_indexed=rules.index_genes_bed.output.bed_gz,
        genes_indexed_tbi=rules.index_genes_bed.output.tbi,
    output:
        plots=OUTPUT.select_track_plots(DataScalingTechnique.PER_SAMPLE),
        template=OUTPUT_DIR + "/track_plots/deeptools/" + DataScalingTechnique.PER_SAMPLE.value + "/template.toml",
    params:
        assay=CONFIG.assay.value,
        peak_files=expand("{peak}.gz", peak=OUTPUT.peak_files),
        genes=OUTPUT_DIR + "/resources/genes_indexed.bed.gz" if CONFIG.assay_config.plot_with_plotnado and CONFIG.genome.genes else None,
        regions=CONFIG.assay_config.plotting.coordinates if CONFIG.assay_config.plotting else None,
        plotting_format=CONFIG.assay_config.plotting.file_format if CONFIG.assay_config.plotting else None,
        outdir=lambda wildcards, output: str(os.path.dirname(output.template)),
    resources:
        mem="1.5GB",
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
    container: "library://asmith151/plotnado/plotnado:latest"
    log: OUTPUT_DIR + "/logs/visualise/plotnado_deeptools_per_sample.log",
    benchmark: OUTPUT_DIR + "/.benchmark/visualise/plotnado_deeptools_per_sample.tsv",
    message: "Generating deeptools per_sample genome browser visualisations with Plotnado"
    script:
        "../../scripts/run_plotnado.py"


# ─── deeptools · individual · per_group ───────────────────────────────────────────
use rule plotnado_deeptools as plotnado_deeptools_per_group with:
    input:
        data=OUTPUT.select_bigwig_subtype(
            method=PileupMethod.DEEPTOOLS,
            scale=DataScalingTechnique.PER_GROUP,
            scaling_method=_DEFAULT_SCALING_METHOD,
            ip_only=True,
        ),
        peaks_indexed=expand("{peak}.gz", peak=OUTPUT.peak_files),
        peaks_tbi=expand("{peak}.gz.tbi", peak=OUTPUT.peak_files),
        genes_indexed=rules.index_genes_bed.output.bed_gz,
        genes_indexed_tbi=rules.index_genes_bed.output.tbi,
    output:
        plots=OUTPUT.select_track_plots(DataScalingTechnique.PER_GROUP, scaling_method=_DEFAULT_SCALING_METHOD),
        template=OUTPUT_DIR + f"/track_plots/deeptools/" + DataScalingTechnique.PER_GROUP.value + f"/{_DEFAULT_SCALING_METHOD}/template.toml",
    log: OUTPUT_DIR + "/logs/visualise/plotnado_deeptools_per_group.log",
    benchmark: OUTPUT_DIR + "/.benchmark/visualise/plotnado_deeptools_per_group.tsv",
    message: "Generating deeptools per-group-scaled genome browser visualisations with Plotnado"


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
        plots=OUTPUT.select_track_plots(DataScalingTechnique.SPIKEIN, spikein_method=SpikeInMethod.ORLANDO.value),
        template=OUTPUT_DIR + "/track_plots/deeptools/spikein/orlando/template.toml",
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
        plots=OUTPUT.select_track_plots(DataScalingTechnique.SPIKEIN, spikein_method=SpikeInMethod.WITH_INPUT.value),
        template=OUTPUT_DIR + "/track_plots/deeptools/spikein/with_input/template.toml",
    log: OUTPUT_DIR + "/logs/visualise/plotnado_deeptools_spikein_withinput.log",
    benchmark: OUTPUT_DIR + "/.benchmark/visualise/plotnado_deeptools_spikein_withinput.tsv",
    message: "Generating deeptools with_input spike-in genome browser visualisations with Plotnado"


# ─── deeptools · merged · per_sample ───────────────────────────────────────────
use rule plotnado_deeptools as plotnado_deeptools_merged with:
    input:
        data=OUTPUT.select_bigwig_subtype(
            method=PileupMethod.DEEPTOOLS,
            scale=DataScalingTechnique.PER_SAMPLE,
            is_merged=True,
        ),
        peaks_indexed=expand("{peak}.gz", peak=OUTPUT.peak_files),
        peaks_tbi=expand("{peak}.gz.tbi", peak=OUTPUT.peak_files),
        genes_indexed=rules.index_genes_bed.output.bed_gz,
        genes_indexed_tbi=rules.index_genes_bed.output.tbi,
    output:
        plots=OUTPUT.select_track_plots(DataScalingTechnique.PER_SAMPLE, is_merged=True),
        template=OUTPUT_DIR + "/track_plots/merged/deeptools/" + DataScalingTechnique.PER_SAMPLE.value + "/template.toml",
    log: OUTPUT_DIR + "/logs/visualise/plotnado_merged_deeptools_per_sample.log",
    benchmark: OUTPUT_DIR + "/.benchmark/visualise/plotnado_merged_deeptools_per_sample.tsv",
    message: "Generating deeptools merged per_sample genome browser visualisations with Plotnado"


# ─── deeptools · merged · per_group ───────────────────────────────────────────────
use rule plotnado_deeptools as plotnado_deeptools_merged_per_group with:
    input:
        data=OUTPUT.select_bigwig_subtype(
            method=PileupMethod.DEEPTOOLS,
            scale=DataScalingTechnique.PER_GROUP,
            scaling_method=_DEFAULT_SCALING_METHOD,
            is_merged=True,
        ),
        peaks_indexed=expand("{peak}.gz", peak=OUTPUT.peak_files),
        peaks_tbi=expand("{peak}.gz.tbi", peak=OUTPUT.peak_files),
        genes_indexed=rules.index_genes_bed.output.bed_gz,
        genes_indexed_tbi=rules.index_genes_bed.output.tbi,
    output:
        plots=OUTPUT.select_track_plots(DataScalingTechnique.PER_GROUP, is_merged=True, scaling_method=_DEFAULT_SCALING_METHOD),
        template=OUTPUT_DIR + f"/track_plots/merged/deeptools/" + DataScalingTechnique.PER_GROUP.value + f"/{_DEFAULT_SCALING_METHOD}/template.toml",
    log: OUTPUT_DIR + "/logs/visualise/plotnado_merged_deeptools_per_group.log",
    benchmark: OUTPUT_DIR + "/.benchmark/visualise/plotnado_merged_deeptools_per_group.tsv",
    message: "Generating deeptools merged per-group-scaled genome browser visualisations with Plotnado"


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
        plots=OUTPUT.select_track_plots(DataScalingTechnique.SPIKEIN, is_merged=True, spikein_method=SpikeInMethod.ORLANDO.value),
        template=OUTPUT_DIR + "/track_plots/merged/deeptools/spikein/orlando/template.toml",
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
        plots=OUTPUT.select_track_plots(DataScalingTechnique.SPIKEIN, is_merged=True, spikein_method=SpikeInMethod.WITH_INPUT.value),
        template=OUTPUT_DIR + "/track_plots/merged/deeptools/spikein/with_input/template.toml",
    log: OUTPUT_DIR + "/logs/visualise/plotnado_merged_deeptools_spikein_withinput.log",
    benchmark: OUTPUT_DIR + "/.benchmark/visualise/plotnado_merged_deeptools_spikein_withinput.tsv",
    message: "Generating deeptools merged with_input spike-in genome browser visualisations with Plotnado"


# ─── bamnado · individual · per_sample ─────────────────────────────────────────
use rule plotnado_deeptools as plotnado_bamnado with:
    input:
        data=OUTPUT.select_bigwig_subtype(
            method=PileupMethod.BAMNADO,
            scale=DataScalingTechnique.PER_SAMPLE,
            ip_only=True,
        ),
        peaks_indexed=expand("{peak}.gz", peak=OUTPUT.peak_files),
        peaks_tbi=expand("{peak}.gz.tbi", peak=OUTPUT.peak_files),
        genes_indexed=rules.index_genes_bed.output.bed_gz,
        genes_indexed_tbi=rules.index_genes_bed.output.tbi,
    output:
        plots=OUTPUT.select_track_plots(DataScalingTechnique.PER_SAMPLE, method=PileupMethod.BAMNADO),
        template=OUTPUT_DIR + "/track_plots/bamnado/" + DataScalingTechnique.PER_SAMPLE.value + "/template.toml",
    log: OUTPUT_DIR + "/logs/visualise/plotnado_bamnado_per_sample.log",
    benchmark: OUTPUT_DIR + "/.benchmark/visualise/plotnado_bamnado_per_sample.tsv",
    message: "Generating bamnado per_sample genome browser visualisations with Plotnado"


# ─── bamnado · merged · per_sample ─────────────────────────────────────────────
use rule plotnado_deeptools as plotnado_bamnado_merged with:
    input:
        data=OUTPUT.select_bigwig_subtype(
            method=PileupMethod.BAMNADO,
            scale=DataScalingTechnique.PER_SAMPLE,
            is_merged=True,
        ),
        peaks_indexed=expand("{peak}.gz", peak=OUTPUT.peak_files),
        peaks_tbi=expand("{peak}.gz.tbi", peak=OUTPUT.peak_files),
        genes_indexed=rules.index_genes_bed.output.bed_gz,
        genes_indexed_tbi=rules.index_genes_bed.output.tbi,
    output:
        plots=OUTPUT.select_track_plots(DataScalingTechnique.PER_SAMPLE, method=PileupMethod.BAMNADO, is_merged=True),
        template=OUTPUT_DIR + "/track_plots/merged/bamnado/" + DataScalingTechnique.PER_SAMPLE.value + "/template.toml",
    log: OUTPUT_DIR + "/logs/visualise/plotnado_merged_bamnado_per_sample.log",
    benchmark: OUTPUT_DIR + "/.benchmark/visualise/plotnado_merged_bamnado_per_sample.tsv",
    message: "Generating bamnado merged per_sample genome browser visualisations with Plotnado"


# ─── bamnado · merged · per_group ─────────────────────────────────────────────────
use rule plotnado_deeptools as plotnado_bamnado_merged_per_group with:
    input:
        data=OUTPUT.select_bigwig_subtype(
            method=PileupMethod.BAMNADO,
            scale=DataScalingTechnique.PER_GROUP,
            scaling_method=_DEFAULT_SCALING_METHOD,
            is_merged=True,
        ),
        peaks_indexed=expand("{peak}.gz", peak=OUTPUT.peak_files),
        peaks_tbi=expand("{peak}.gz.tbi", peak=OUTPUT.peak_files),
        genes_indexed=rules.index_genes_bed.output.bed_gz,
        genes_indexed_tbi=rules.index_genes_bed.output.tbi,
    output:
        plots=OUTPUT.select_track_plots(DataScalingTechnique.PER_GROUP, method=PileupMethod.BAMNADO, is_merged=True, scaling_method=_DEFAULT_SCALING_METHOD),
        template=OUTPUT_DIR + f"/track_plots/merged/bamnado/" + DataScalingTechnique.PER_GROUP.value + f"/{_DEFAULT_SCALING_METHOD}/template.toml",
    log: OUTPUT_DIR + "/logs/visualise/plotnado_merged_bamnado_per_group.log",
    benchmark: OUTPUT_DIR + "/.benchmark/visualise/plotnado_merged_bamnado_per_group.tsv",
    message: "Generating bamnado merged per-group-scaled genome browser visualisations with Plotnado"


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
        plots=OUTPUT.select_track_plots(DataScalingTechnique.SPIKEIN, method=PileupMethod.BAMNADO, is_merged=True, spikein_method=SpikeInMethod.ORLANDO.value),
        template=OUTPUT_DIR + "/track_plots/merged/bamnado/spikein/orlando/template.toml",
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
        plots=OUTPUT.select_track_plots(DataScalingTechnique.SPIKEIN, method=PileupMethod.BAMNADO, is_merged=True, spikein_method=SpikeInMethod.WITH_INPUT.value),
        template=OUTPUT_DIR + "/track_plots/merged/bamnado/spikein/with_input/template.toml",
    log: OUTPUT_DIR + "/logs/visualise/plotnado_merged_bamnado_spikein_withinput.log",
    benchmark: OUTPUT_DIR + "/.benchmark/visualise/plotnado_merged_bamnado_spikein_withinput.tsv",
    message: "Generating bamnado merged with_input spike-in genome browser visualisations with Plotnado"


# ─── homer · individual · per_sample ───────────────────────────────────────────
use rule plotnado_deeptools as plotnado_homer with:
    input:
        data=OUTPUT.select_bigwig_subtype(
            method=PileupMethod.HOMER,
            scale=DataScalingTechnique.PER_SAMPLE,
            ip_only=True,
        ),
        peaks_indexed=expand("{peak}.gz", peak=OUTPUT.peak_files),
        peaks_tbi=expand("{peak}.gz.tbi", peak=OUTPUT.peak_files),
        genes_indexed=rules.index_genes_bed.output.bed_gz,
        genes_indexed_tbi=rules.index_genes_bed.output.tbi,
    output:
        plots=OUTPUT.select_track_plots(DataScalingTechnique.PER_SAMPLE, method=PileupMethod.HOMER),
        template=OUTPUT_DIR + "/track_plots/homer/" + DataScalingTechnique.PER_SAMPLE.value + "/template.toml",
    log: OUTPUT_DIR + "/logs/visualise/plotnado_homer_per_sample.log",
    benchmark: OUTPUT_DIR + "/.benchmark/visualise/plotnado_homer_per_sample.tsv",
    message: "Generating homer per_sample genome browser visualisations with Plotnado"


# ─── homer · merged · per_sample ───────────────────────────────────────────────
use rule plotnado_deeptools as plotnado_homer_merged with:
    input:
        data=OUTPUT.select_bigwig_subtype(
            method=PileupMethod.HOMER,
            scale=DataScalingTechnique.PER_SAMPLE,
            is_merged=True,
        ),
        peaks_indexed=expand("{peak}.gz", peak=OUTPUT.peak_files),
        peaks_tbi=expand("{peak}.gz.tbi", peak=OUTPUT.peak_files),
        genes_indexed=rules.index_genes_bed.output.bed_gz,
        genes_indexed_tbi=rules.index_genes_bed.output.tbi,
    output:
        plots=OUTPUT.select_track_plots(DataScalingTechnique.PER_SAMPLE, method=PileupMethod.HOMER, is_merged=True),
        template=OUTPUT_DIR + "/track_plots/merged/homer/" + DataScalingTechnique.PER_SAMPLE.value + "/template.toml",
    log: OUTPUT_DIR + "/logs/visualise/plotnado_merged_homer_per_sample.log",
    benchmark: OUTPUT_DIR + "/.benchmark/visualise/plotnado_merged_homer_per_sample.tsv",
    message: "Generating homer merged per_sample genome browser visualisations with Plotnado"


# ─── methyldackel · individual · per_sample (METH assay only) ──────────────────
if ASSAY == Assay.METH:
    _meth_method = CONFIG.assay_config.methylation.method.value  # "taps" or "wgbs"

    rule plotnado_methyldackel:
        input:
            data=OUTPUT.select_files(".bigWig", include=f"bigwigs/{_meth_method}/", exclude="geo_submission"),
            genes_indexed=rules.index_genes_bed.output.bed_gz,
            genes_indexed_tbi=rules.index_genes_bed.output.tbi,
        output:
            plots=OUTPUT.select_track_plots(DataScalingTechnique.PER_SAMPLE, method=PileupMethod.METHYLDACKEL),
            template=OUTPUT_DIR + "/track_plots/methyldackel/" + DataScalingTechnique.PER_SAMPLE.value + "/template.toml",
        params:
            assay=CONFIG.assay.value,
            peak_files=[],
            genes=OUTPUT_DIR + "/resources/genes_indexed.bed.gz" if CONFIG.assay_config.plot_with_plotnado and CONFIG.genome.genes else None,
            regions=CONFIG.assay_config.plotting.coordinates if CONFIG.assay_config.plotting else None,
            plotting_format=CONFIG.assay_config.plotting.file_format if CONFIG.assay_config.plotting else None,
            outdir=lambda wildcards, output: str(os.path.dirname(output.template)),
        resources:
            mem="1.5GB",
            runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
        container: "library://asmith151/plotnado/plotnado:latest"
        log: OUTPUT_DIR + "/logs/visualise/plotnado_methyldackel_per_sample.log"
        benchmark: OUTPUT_DIR + "/.benchmark/visualise/plotnado_methyldackel_per_sample.tsv"
        message: "Generating methyldackel per_sample genome browser visualisations with Plotnado"
        script:
            "../../scripts/run_plotnado.py"
