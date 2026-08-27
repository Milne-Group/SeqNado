from seqnado import Assay, PileupMethod, DataScalingTechnique, GroupScalingMethod, SpikeInMethod
from seqnado.workflow.helpers.common import define_memory_requested, define_time_requested

_spikein_cfg = CONFIG.assay_config.spikein
_spikein_methods = _spikein_cfg.method if _spikein_cfg else []
_has_spikein_orlando = SpikeInMethod.ORLANDO in _spikein_methods
_has_spikein_withinput = SpikeInMethod.WITH_INPUT in _spikein_methods

# Heatmaps for the per-group technique are generated for one bamnado scaling
# method only — additional configured group_scaling_methods still produce
# bigwigs, just not dedicated heatmaps. OUTPUT picks which one, so the rules and
# the expected-file list cannot drift apart.
_DEFAULT_SCALING_METHOD = OUTPUT.plot_scaling_method.value

# Heatmap and metaplot generation from bigWig files.
# One set of rules per (pileup method) × (scale) × (merged) combination.
# Compatible combinations:
#   deeptools    : individual → per_sample/per_group/spikein_orlando/spikein_withinput ; merged → per_sample/per_group/spikein_orlando/spikein_withinput
#   bamnado      : individual → per_sample/spikein_orlando/spikein_withinput      ; merged → per_sample/per_group/spikein_orlando/spikein_withinput
#   (per_group is genomics-only — excluded entirely for RNA-seq)
#   homer        : individual → per_sample              ; merged → per_sample
#   methyldackel : individual → per_sample (METH assay only, reads from bigwigs/taps or bigwigs/wgbs)


# ─── methyldackel · individual · per_sample (METH assay only) ──────────────────
if ASSAY == Assay.METH:
    _meth_method = CONFIG.assay_config.methylation.method.value  # "taps" or "wgbs"

    _deeptools_matrix_options = str(CONFIG.third_party_tools.deeptools.compute_matrix.command_line_arguments) if CONFIG.third_party_tools.deeptools else ""
    _deeptools_matrix_threads = CONFIG.third_party_tools.deeptools.compute_matrix.threads if CONFIG.third_party_tools.deeptools else 4
    _deeptools_heatmap_options = str(CONFIG.third_party_tools.deeptools.plot_heatmap.command_line_arguments) if CONFIG.third_party_tools.deeptools else ""

    rule heatmap_methyldackel_per_sample_matrix:
        input:
            bigwigs=OUTPUT.select_files(".bigWig", include=f"bigwigs/{_meth_method}/", exclude="geo_submission"),
        output:
            matrix=temp(OUTPUT.select_heatmap_matrix(
                DataScalingTechnique.PER_SAMPLE,
                method=PileupMethod.METHYLDACKEL,
            )),
        params:
            gtf=CONFIG.genome.gtf,
            options=_deeptools_matrix_options,
        threads: _deeptools_matrix_threads
        resources:
            runtime=lambda wildcards, attempt: f"{1 * 2**attempt}h",
            mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
        log: OUTPUT_DIR + "/logs/heatmap/methyldackel/" + DataScalingTechnique.PER_SAMPLE.value + "/matrix.log"
        benchmark: OUTPUT_DIR + "/.benchmark/heatmap/methyldackel/" + DataScalingTechnique.PER_SAMPLE.value + "/matrix.tsv"
        message: "Computing methylation heatmap matrix from bigWig files"
        shell: """
        computeMatrix scale-regions -p {threads} {params.options} --smartLabels --missingDataAsZero -S {input.bigwigs} -R {params.gtf} -o {output.matrix} >> {log} 2>&1
        """

    rule heatmap_methyldackel_per_sample_plot:
        input:
            matrix=rules.heatmap_methyldackel_per_sample_matrix.output.matrix,
        output:
            heatmap=OUTPUT.select_heatmap_plot(
                DataScalingTechnique.PER_SAMPLE,
                method=PileupMethod.METHYLDACKEL,
            ),
        params:
            options=_deeptools_heatmap_options,
        resources:
            mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
        log: OUTPUT_DIR + "/logs/heatmap/methyldackel/" + DataScalingTechnique.PER_SAMPLE.value + "/heatmap.log"
        benchmark: OUTPUT_DIR + "/.benchmark/heatmap/methyldackel/" + DataScalingTechnique.PER_SAMPLE.value + "/heatmap.tsv"
        message: "Generating methylation heatmap from matrix"
        shell: """
        plotHeatmap -m {input.matrix} -out {output.heatmap} {params.options}
        """

    rule heatmap_methyldackel_per_sample_metaplot:
        input:
            matrix=rules.heatmap_methyldackel_per_sample_matrix.output.matrix,
        output:
            metaplot=OUTPUT.select_heatmap_metaplot(
                DataScalingTechnique.PER_SAMPLE,
                method=PileupMethod.METHYLDACKEL,
            ),
        resources:
            mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
        log: OUTPUT_DIR + "/logs/heatmap/methyldackel/" + DataScalingTechnique.PER_SAMPLE.value + "/metaplot.log"
        benchmark: OUTPUT_DIR + "/.benchmark/heatmap/methyldackel/" + DataScalingTechnique.PER_SAMPLE.value + "/metaplot.tsv"
        message: "Generating methylation metaplot from heatmap matrix"
        shell: """
        plotProfile -m {input.matrix} -out {output.metaplot} --perGroup
        """


# ─── deeptools · individual · per_sample (base rules) ──────────────────────────
rule heatmap_deeptools_per_sample_matrix:
    input:
        bigwigs=OUTPUT.select_bigwig_subtype(
            method=PileupMethod.DEEPTOOLS,
            scale=DataScalingTechnique.PER_SAMPLE,
            ip_only=True,
        ),
    output:
        matrix=temp(OUTPUT.select_heatmap_matrix(
            DataScalingTechnique.PER_SAMPLE
            )),
    params:
        gtf=CONFIG.genome.gtf,
        options=str(CONFIG.third_party_tools.deeptools.compute_matrix.command_line_arguments),
    threads: CONFIG.third_party_tools.deeptools.compute_matrix.threads
    resources:
        runtime=lambda wildcards, attempt: f"{1 * 2**attempt}h",
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/heatmap/deeptools/" + DataScalingTechnique.PER_SAMPLE.value + "/matrix.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/deeptools/" + DataScalingTechnique.PER_SAMPLE.value + "/matrix.tsv",
    message: "Computing deeptools per_sample heatmap matrix from bigWig files"
    shell: """
    computeMatrix scale-regions -p {threads} {params.options} --smartLabels --missingDataAsZero -S {input.bigwigs} -R {params.gtf} -o {output.matrix} >> {log} 2>&1
    """


rule heatmap_deeptools_per_sample_plot:
    input:
        matrix=rules.heatmap_deeptools_per_sample_matrix.output.matrix,
    output:
        heatmap=OUTPUT.select_heatmap_plot(
            DataScalingTechnique.PER_SAMPLE),
    params:
        options=str(CONFIG.third_party_tools.deeptools.plot_heatmap.command_line_arguments),
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/heatmap/deeptools/" + DataScalingTechnique.PER_SAMPLE.value + "/heatmap.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/deeptools/" + DataScalingTechnique.PER_SAMPLE.value + "/heatmap.tsv",
    message: "Generating deeptools per_sample heatmap from matrix"
    shell: """
    plotHeatmap -m {input.matrix} -out {output.heatmap} {params.options}
    """


rule heatmap_deeptools_per_sample_metaplot:
    input:
        matrix=rules.heatmap_deeptools_per_sample_matrix.output.matrix,
    output:
        metaplot=OUTPUT.select_heatmap_metaplot(
            DataScalingTechnique.PER_SAMPLE),
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/heatmap/deeptools/" + DataScalingTechnique.PER_SAMPLE.value + "/metaplot.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/deeptools/" + DataScalingTechnique.PER_SAMPLE.value + "/metaplot.tsv",
    message: "Generating deeptools per_sample metaplot from heatmap matrix"
    shell: """
    plotProfile -m {input.matrix} -out {output.metaplot} --perGroup
    """


# ─── deeptools · individual · per_group ───────────────────────────────────────────
use rule heatmap_deeptools_per_sample_matrix as heatmap_deeptools_per_group_matrix with:
    input:
        bigwigs=OUTPUT.select_bigwig_subtype(
            method=PileupMethod.DEEPTOOLS,
            scale=DataScalingTechnique.PER_GROUP,
            scaling_method=_DEFAULT_SCALING_METHOD,
            ip_only=True,
        ),
    output:
        matrix=temp(OUTPUT.select_heatmap_matrix(
            DataScalingTechnique.PER_GROUP, scaling_method=_DEFAULT_SCALING_METHOD)),
    log: OUTPUT_DIR + "/logs/heatmap/deeptools/" + DataScalingTechnique.PER_GROUP.value + "/matrix.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/deeptools/" + DataScalingTechnique.PER_GROUP.value + "/matrix.tsv",
    message: "Computing deeptools per-group-scaled heatmap matrix from bigWig files"

use rule heatmap_deeptools_per_sample_plot as heatmap_deeptools_per_group_plot with:
    input:
        matrix=OUTPUT.select_heatmap_matrix(
            DataScalingTechnique.PER_GROUP, scaling_method=_DEFAULT_SCALING_METHOD),
    output:
        heatmap=OUTPUT.select_heatmap_plot(
            DataScalingTechnique.PER_GROUP, scaling_method=_DEFAULT_SCALING_METHOD),
    log: OUTPUT_DIR + "/logs/heatmap/deeptools/" + DataScalingTechnique.PER_GROUP.value + "/heatmap.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/deeptools/" + DataScalingTechnique.PER_GROUP.value + "/heatmap.tsv",
    message: "Generating deeptools per-group-scaled heatmap from matrix"

use rule heatmap_deeptools_per_sample_metaplot as heatmap_deeptools_per_group_metaplot with:
    input:
        matrix=OUTPUT.select_heatmap_matrix(
            DataScalingTechnique.PER_GROUP, scaling_method=_DEFAULT_SCALING_METHOD),
    output:
        metaplot=OUTPUT.select_heatmap_metaplot(
            DataScalingTechnique.PER_GROUP, scaling_method=_DEFAULT_SCALING_METHOD),
    log: OUTPUT_DIR + "/logs/heatmap/deeptools/" + DataScalingTechnique.PER_GROUP.value + "/metaplot.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/deeptools/" + DataScalingTechnique.PER_GROUP.value + "/metaplot.tsv",
    message: "Generating deeptools per-group-scaled metaplot from matrix"


# ─── deeptools · individual · spikein · orlando ───────────────────────────────
if _has_spikein_orlando:
    use rule heatmap_deeptools_per_sample_matrix as heatmap_deeptools_spikein_orlando_matrix with:
        input:
            bigwigs=OUTPUT.select_bigwig_subtype(
                method=PileupMethod.DEEPTOOLS,
                scale=DataScalingTechnique.SPIKEIN,
                spikein_method=SpikeInMethod.ORLANDO.value,
                ip_only=True,
            ),
        output:
            matrix=temp(OUTPUT.select_heatmap_matrix(
            DataScalingTechnique.SPIKEIN,
            spikein_method=SpikeInMethod.ORLANDO.value)),
        log: OUTPUT_DIR + "/logs/heatmap/deeptools/spikein/orlando/matrix.log",
        benchmark: OUTPUT_DIR + "/.benchmark/heatmap/deeptools/spikein/orlando/matrix.tsv",
        message: "Computing deeptools orlando spike-in heatmap matrix from bigWig files"

    use rule heatmap_deeptools_per_sample_plot as heatmap_deeptools_spikein_orlando_plot with:
        input:
            matrix=OUTPUT.select_heatmap_matrix(
            DataScalingTechnique.SPIKEIN,
            spikein_method=SpikeInMethod.ORLANDO.value),
        output:
            heatmap=OUTPUT.select_heatmap_plot(
            DataScalingTechnique.SPIKEIN,
            spikein_method=SpikeInMethod.ORLANDO.value),
        log: OUTPUT_DIR + "/logs/heatmap/deeptools/spikein/orlando/heatmap.log",
        benchmark: OUTPUT_DIR + "/.benchmark/heatmap/deeptools/spikein/orlando/heatmap.tsv",
        message: "Generating deeptools orlando spike-in heatmap from matrix"

    use rule heatmap_deeptools_per_sample_metaplot as heatmap_deeptools_spikein_orlando_metaplot with:
        input:
            matrix=OUTPUT.select_heatmap_matrix(
            DataScalingTechnique.SPIKEIN,
            spikein_method=SpikeInMethod.ORLANDO.value),
        output:
            metaplot=OUTPUT.select_heatmap_metaplot(
            DataScalingTechnique.SPIKEIN,
            spikein_method=SpikeInMethod.ORLANDO.value),
        log: OUTPUT_DIR + "/logs/heatmap/deeptools/spikein/orlando/metaplot.log",
        benchmark: OUTPUT_DIR + "/.benchmark/heatmap/deeptools/spikein/orlando/metaplot.tsv",
        message: "Generating deeptools orlando spike-in metaplot from matrix"


# ─── deeptools · individual · spikein · with_input ───────────────────────────
if _has_spikein_withinput:
    use rule heatmap_deeptools_per_sample_matrix as heatmap_deeptools_spikein_withinput_matrix with:
        input:
            bigwigs=OUTPUT.select_bigwig_subtype(
                method=PileupMethod.DEEPTOOLS,
                scale=DataScalingTechnique.SPIKEIN,
                spikein_method=SpikeInMethod.WITH_INPUT.value,
                ip_only=True,
            ),
        output:
            matrix=temp(OUTPUT.select_heatmap_matrix(
            DataScalingTechnique.SPIKEIN,
            spikein_method=SpikeInMethod.WITH_INPUT.value)),
        log: OUTPUT_DIR + "/logs/heatmap/deeptools/spikein/with_input/matrix.log",
        benchmark: OUTPUT_DIR + "/.benchmark/heatmap/deeptools/spikein/with_input/matrix.tsv",
        message: "Computing deeptools with_input spike-in heatmap matrix from bigWig files"

    use rule heatmap_deeptools_per_sample_plot as heatmap_deeptools_spikein_withinput_plot with:
        input:
            matrix=OUTPUT.select_heatmap_matrix(
            DataScalingTechnique.SPIKEIN,
            spikein_method=SpikeInMethod.WITH_INPUT.value),
        output:
            heatmap=OUTPUT.select_heatmap_plot(
            DataScalingTechnique.SPIKEIN,
            spikein_method=SpikeInMethod.WITH_INPUT.value),
        log: OUTPUT_DIR + "/logs/heatmap/deeptools/spikein/with_input/heatmap.log",
        benchmark: OUTPUT_DIR + "/.benchmark/heatmap/deeptools/spikein/with_input/heatmap.tsv",
        message: "Generating deeptools with_input spike-in heatmap from matrix"

    use rule heatmap_deeptools_per_sample_metaplot as heatmap_deeptools_spikein_withinput_metaplot with:
        input:
            matrix=OUTPUT.select_heatmap_matrix(
            DataScalingTechnique.SPIKEIN,
            spikein_method=SpikeInMethod.WITH_INPUT.value),
        output:
            metaplot=OUTPUT.select_heatmap_metaplot(
            DataScalingTechnique.SPIKEIN,
            spikein_method=SpikeInMethod.WITH_INPUT.value),
        log: OUTPUT_DIR + "/logs/heatmap/deeptools/spikein/with_input/metaplot.log",
        benchmark: OUTPUT_DIR + "/.benchmark/heatmap/deeptools/spikein/with_input/metaplot.tsv",
        message: "Generating deeptools with_input spike-in metaplot from matrix"


# ─── deeptools · merged · per_sample ───────────────────────────────────────────
use rule heatmap_deeptools_per_sample_matrix as heatmap_deeptools_merged_matrix with:
    input:
        bigwigs=OUTPUT.select_bigwig_subtype(
            method=PileupMethod.DEEPTOOLS,
            scale=DataScalingTechnique.PER_SAMPLE,
            is_merged=True,
        ),
    output:
        matrix=temp(OUTPUT.select_heatmap_matrix(
            DataScalingTechnique.PER_SAMPLE,
            is_merged=True)),
    log: OUTPUT_DIR + "/logs/heatmap/merged/deeptools/" + DataScalingTechnique.PER_SAMPLE.value + "/matrix.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/merged/deeptools/" + DataScalingTechnique.PER_SAMPLE.value + "/matrix.tsv",
    message: "Computing deeptools merged per_sample heatmap matrix from bigWig files"

use rule heatmap_deeptools_per_sample_plot as heatmap_deeptools_merged_plot with:
    input:
        matrix=OUTPUT.select_heatmap_matrix(
            DataScalingTechnique.PER_SAMPLE,
            is_merged=True),
    output:
        heatmap=OUTPUT.select_heatmap_plot(
            DataScalingTechnique.PER_SAMPLE,
            is_merged=True),
    log: OUTPUT_DIR + "/logs/heatmap/merged/deeptools/" + DataScalingTechnique.PER_SAMPLE.value + "/heatmap.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/merged/deeptools/" + DataScalingTechnique.PER_SAMPLE.value + "/heatmap.tsv",
    message: "Generating deeptools merged per_sample heatmap from matrix"

use rule heatmap_deeptools_per_sample_metaplot as heatmap_deeptools_merged_metaplot with:
    input:
        matrix=OUTPUT.select_heatmap_matrix(
            DataScalingTechnique.PER_SAMPLE,
            is_merged=True),
    output:
        metaplot=OUTPUT.select_heatmap_metaplot(
            DataScalingTechnique.PER_SAMPLE,
            is_merged=True),
    log: OUTPUT_DIR + "/logs/heatmap/merged/deeptools/" + DataScalingTechnique.PER_SAMPLE.value + "/metaplot.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/merged/deeptools/" + DataScalingTechnique.PER_SAMPLE.value + "/metaplot.tsv",
    message: "Generating deeptools merged per_sample metaplot from matrix"


# ─── deeptools · merged · per_group ───────────────────────────────────────────────
use rule heatmap_deeptools_per_sample_matrix as heatmap_deeptools_merged_per_group_matrix with:
    input:
        bigwigs=OUTPUT.select_bigwig_subtype(
            method=PileupMethod.DEEPTOOLS,
            scale=DataScalingTechnique.PER_GROUP,
            scaling_method=_DEFAULT_SCALING_METHOD,
            is_merged=True,
        ),
    output:
        matrix=temp(OUTPUT.select_heatmap_matrix(
            DataScalingTechnique.PER_GROUP,
            is_merged=True, scaling_method=_DEFAULT_SCALING_METHOD)),
    log: OUTPUT_DIR + "/logs/heatmap/merged/deeptools/" + DataScalingTechnique.PER_GROUP.value + "/matrix.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/merged/deeptools/" + DataScalingTechnique.PER_GROUP.value + "/matrix.tsv",
    message: "Computing deeptools merged per-group-scaled heatmap matrix from bigWig files"

use rule heatmap_deeptools_per_sample_plot as heatmap_deeptools_merged_per_group_plot with:
    input:
        matrix=OUTPUT.select_heatmap_matrix(
            DataScalingTechnique.PER_GROUP,
            is_merged=True, scaling_method=_DEFAULT_SCALING_METHOD),
    output:
        heatmap=OUTPUT.select_heatmap_plot(
            DataScalingTechnique.PER_GROUP,
            is_merged=True, scaling_method=_DEFAULT_SCALING_METHOD),
    log: OUTPUT_DIR + "/logs/heatmap/merged/deeptools/" + DataScalingTechnique.PER_GROUP.value + "/heatmap.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/merged/deeptools/" + DataScalingTechnique.PER_GROUP.value + "/heatmap.tsv",
    message: "Generating deeptools merged per-group-scaled heatmap from matrix"

use rule heatmap_deeptools_per_sample_metaplot as heatmap_deeptools_merged_per_group_metaplot with:
    input:
        matrix=OUTPUT.select_heatmap_matrix(
            DataScalingTechnique.PER_GROUP,
            is_merged=True, scaling_method=_DEFAULT_SCALING_METHOD),
    output:
        metaplot=OUTPUT.select_heatmap_metaplot(
            DataScalingTechnique.PER_GROUP,
            is_merged=True, scaling_method=_DEFAULT_SCALING_METHOD),
    log: OUTPUT_DIR + "/logs/heatmap/merged/deeptools/" + DataScalingTechnique.PER_GROUP.value + "/metaplot.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/merged/deeptools/" + DataScalingTechnique.PER_GROUP.value + "/metaplot.tsv",
    message: "Generating deeptools merged per-group-scaled metaplot from matrix"


# ─── deeptools · merged · spikein · orlando ──────────────────────────────────
if _has_spikein_orlando:
    use rule heatmap_deeptools_per_sample_matrix as heatmap_deeptools_merged_spikein_orlando_matrix with:
        input:
            bigwigs=OUTPUT.select_bigwig_subtype(
                method=PileupMethod.DEEPTOOLS,
                scale=DataScalingTechnique.SPIKEIN,
                spikein_method=SpikeInMethod.ORLANDO.value,
                is_merged=True,
            ),
        output:
            matrix=temp(OUTPUT.select_heatmap_matrix(
            DataScalingTechnique.SPIKEIN,
            is_merged=True,
            spikein_method=SpikeInMethod.ORLANDO.value)),
        log: OUTPUT_DIR + "/logs/heatmap/merged/deeptools/spikein/orlando/matrix.log",
        benchmark: OUTPUT_DIR + "/.benchmark/heatmap/merged/deeptools/spikein/orlando/matrix.tsv",
        message: "Computing deeptools merged orlando spike-in heatmap matrix from bigWig files"

    use rule heatmap_deeptools_per_sample_plot as heatmap_deeptools_merged_spikein_orlando_plot with:
        input:
            matrix=OUTPUT.select_heatmap_matrix(
            DataScalingTechnique.SPIKEIN,
            is_merged=True,
            spikein_method=SpikeInMethod.ORLANDO.value),
        output:
            heatmap=OUTPUT.select_heatmap_plot(
            DataScalingTechnique.SPIKEIN,
            is_merged=True,
            spikein_method=SpikeInMethod.ORLANDO.value),
        log: OUTPUT_DIR + "/logs/heatmap/merged/deeptools/spikein/orlando/heatmap.log",
        benchmark: OUTPUT_DIR + "/.benchmark/heatmap/merged/deeptools/spikein/orlando/heatmap.tsv",
        message: "Generating deeptools merged orlando spike-in heatmap from matrix"

    use rule heatmap_deeptools_per_sample_metaplot as heatmap_deeptools_merged_spikein_orlando_metaplot with:
        input:
            matrix=OUTPUT.select_heatmap_matrix(
            DataScalingTechnique.SPIKEIN,
            is_merged=True,
            spikein_method=SpikeInMethod.ORLANDO.value),
        output:
            metaplot=OUTPUT.select_heatmap_metaplot(
            DataScalingTechnique.SPIKEIN,
            is_merged=True,
            spikein_method=SpikeInMethod.ORLANDO.value),
        log: OUTPUT_DIR + "/logs/heatmap/merged/deeptools/spikein/orlando/metaplot.log",
        benchmark: OUTPUT_DIR + "/.benchmark/heatmap/merged/deeptools/spikein/orlando/metaplot.tsv",
        message: "Generating deeptools merged orlando spike-in metaplot from matrix"


# ─── deeptools · merged · spikein · with_input ───────────────────────────────
if _has_spikein_withinput:
    use rule heatmap_deeptools_per_sample_matrix as heatmap_deeptools_merged_spikein_withinput_matrix with:
        input:
            bigwigs=OUTPUT.select_bigwig_subtype(
                method=PileupMethod.DEEPTOOLS,
                scale=DataScalingTechnique.SPIKEIN,
                spikein_method=SpikeInMethod.WITH_INPUT.value,
                is_merged=True,
            ),
        output:
            matrix=temp(OUTPUT.select_heatmap_matrix(
            DataScalingTechnique.SPIKEIN,
            is_merged=True,
            spikein_method=SpikeInMethod.WITH_INPUT.value)),
        log: OUTPUT_DIR + "/logs/heatmap/merged/deeptools/spikein/with_input/matrix.log",
        benchmark: OUTPUT_DIR + "/.benchmark/heatmap/merged/deeptools/spikein/with_input/matrix.tsv",
        message: "Computing deeptools merged with_input spike-in heatmap matrix from bigWig files"

    use rule heatmap_deeptools_per_sample_plot as heatmap_deeptools_merged_spikein_withinput_plot with:
        input:
            matrix=OUTPUT.select_heatmap_matrix(
            DataScalingTechnique.SPIKEIN,
            is_merged=True,
            spikein_method=SpikeInMethod.WITH_INPUT.value),
        output:
            heatmap=OUTPUT.select_heatmap_plot(
            DataScalingTechnique.SPIKEIN,
            is_merged=True,
            spikein_method=SpikeInMethod.WITH_INPUT.value),
        log: OUTPUT_DIR + "/logs/heatmap/merged/deeptools/spikein/with_input/heatmap.log",
        benchmark: OUTPUT_DIR + "/.benchmark/heatmap/merged/deeptools/spikein/with_input/heatmap.tsv",
        message: "Generating deeptools merged with_input spike-in heatmap from matrix"

    use rule heatmap_deeptools_per_sample_metaplot as heatmap_deeptools_merged_spikein_withinput_metaplot with:
        input:
            matrix=OUTPUT.select_heatmap_matrix(
            DataScalingTechnique.SPIKEIN,
            is_merged=True,
            spikein_method=SpikeInMethod.WITH_INPUT.value),
        output:
            metaplot=OUTPUT.select_heatmap_metaplot(
            DataScalingTechnique.SPIKEIN,
            is_merged=True,
            spikein_method=SpikeInMethod.WITH_INPUT.value),
        log: OUTPUT_DIR + "/logs/heatmap/merged/deeptools/spikein/with_input/metaplot.log",
        benchmark: OUTPUT_DIR + "/.benchmark/heatmap/merged/deeptools/spikein/with_input/metaplot.tsv",
        message: "Generating deeptools merged with_input spike-in metaplot from matrix"


# ─── bamnado · individual · per_sample ─────────────────────────────────────────
use rule heatmap_deeptools_per_sample_matrix as heatmap_bamnado_per_sample_matrix with:
    input:
        bigwigs=OUTPUT.select_bigwig_subtype(
            method=PileupMethod.BAMNADO,
            scale=DataScalingTechnique.PER_SAMPLE,
            ip_only=True,
        ),
    output:
        matrix=temp(OUTPUT.select_heatmap_matrix(
            DataScalingTechnique.PER_SAMPLE, 
            method=PileupMethod.BAMNADO)),
    log: OUTPUT_DIR + "/logs/heatmap/bamnado/" + DataScalingTechnique.PER_SAMPLE.value + "/matrix.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/bamnado/" + DataScalingTechnique.PER_SAMPLE.value + "/matrix.tsv",
    message: "Computing bamnado per_sample heatmap matrix from bigWig files"

use rule heatmap_deeptools_per_sample_plot as heatmap_bamnado_per_sample_plot with:
    input:
        matrix=OUTPUT.select_heatmap_matrix(
            DataScalingTechnique.PER_SAMPLE, 
            method=PileupMethod.BAMNADO),
    output:
        heatmap=OUTPUT.select_heatmap_plot(
            DataScalingTechnique.PER_SAMPLE, 
            method=PileupMethod.BAMNADO),
    log: OUTPUT_DIR + "/logs/heatmap/bamnado/" + DataScalingTechnique.PER_SAMPLE.value + "/heatmap.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/bamnado/" + DataScalingTechnique.PER_SAMPLE.value + "/heatmap.tsv",
    message: "Generating bamnado per_sample heatmap from matrix"

use rule heatmap_deeptools_per_sample_metaplot as heatmap_bamnado_per_sample_metaplot with:
    input:
        matrix=OUTPUT.select_heatmap_matrix(
            DataScalingTechnique.PER_SAMPLE, 
            method=PileupMethod.BAMNADO),
    output:
        metaplot=OUTPUT.select_heatmap_metaplot(
            DataScalingTechnique.PER_SAMPLE, 
            method=PileupMethod.BAMNADO),
    log: OUTPUT_DIR + "/logs/heatmap/bamnado/" + DataScalingTechnique.PER_SAMPLE.value + "/metaplot.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/bamnado/" + DataScalingTechnique.PER_SAMPLE.value + "/metaplot.tsv",
    message: "Generating bamnado per_sample metaplot from matrix"


# ─── bamnado · merged · per_sample ─────────────────────────────────────────────
use rule heatmap_deeptools_per_sample_matrix as heatmap_bamnado_merged_per_sample_matrix with:
    input:
        bigwigs=OUTPUT.select_bigwig_subtype(
            method=PileupMethod.BAMNADO,
            scale=DataScalingTechnique.PER_SAMPLE,
            is_merged=True,
        ),
    output:
        matrix=temp(OUTPUT.select_heatmap_matrix(
            DataScalingTechnique.PER_SAMPLE, 
            method=PileupMethod.BAMNADO,
            is_merged=True)),
    log: OUTPUT_DIR + "/logs/heatmap/merged/bamnado/" + DataScalingTechnique.PER_SAMPLE.value + "/matrix.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/merged/bamnado/" + DataScalingTechnique.PER_SAMPLE.value + "/matrix.tsv",
    message: "Computing bamnado merged per_sample heatmap matrix from bigWig files"

use rule heatmap_deeptools_per_sample_plot as heatmap_bamnado_merged_per_sample_plot with:
    input:
        matrix=OUTPUT.select_heatmap_matrix(
            DataScalingTechnique.PER_SAMPLE, 
            method=PileupMethod.BAMNADO,
            is_merged=True),
    output:
        heatmap=OUTPUT.select_heatmap_plot(
            DataScalingTechnique.PER_SAMPLE, 
            method=PileupMethod.BAMNADO,
            is_merged=True),
    log: OUTPUT_DIR + "/logs/heatmap/merged/bamnado/" + DataScalingTechnique.PER_SAMPLE.value + "/heatmap.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/merged/bamnado/" + DataScalingTechnique.PER_SAMPLE.value + "/heatmap.tsv",
    message: "Generating bamnado merged per_sample heatmap from matrix"

use rule heatmap_deeptools_per_sample_metaplot as heatmap_bamnado_merged_per_sample_metaplot with:
    input:
        matrix=OUTPUT.select_heatmap_matrix(
            DataScalingTechnique.PER_SAMPLE, 
            method=PileupMethod.BAMNADO,
            is_merged=True),
    output:
        metaplot=OUTPUT.select_heatmap_metaplot(
            DataScalingTechnique.PER_SAMPLE, 
            method=PileupMethod.BAMNADO,
            is_merged=True),
    log: OUTPUT_DIR + "/logs/heatmap/merged/bamnado/" + DataScalingTechnique.PER_SAMPLE.value + "/metaplot.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/merged/bamnado/" + DataScalingTechnique.PER_SAMPLE.value + "/metaplot.tsv",
    message: "Generating bamnado merged per_sample metaplot from matrix"


# ─── bamnado · merged · per_group ─────────────────────────────────────────────────
use rule heatmap_deeptools_per_sample_matrix as heatmap_bamnado_merged_per_group_matrix with:
    input:
        bigwigs=OUTPUT.select_bigwig_subtype(
            method=PileupMethod.BAMNADO,
            scale=DataScalingTechnique.PER_GROUP,
            scaling_method=_DEFAULT_SCALING_METHOD,
            is_merged=True,
        ),
    output:
        matrix=temp(OUTPUT.select_heatmap_matrix(
            DataScalingTechnique.PER_GROUP,
            method=PileupMethod.BAMNADO,
            is_merged=True, scaling_method=_DEFAULT_SCALING_METHOD)),
    log: OUTPUT_DIR + "/logs/heatmap/merged/bamnado/" + DataScalingTechnique.PER_GROUP.value + "/matrix.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/merged/bamnado/" + DataScalingTechnique.PER_GROUP.value + "/matrix.tsv",
    message: "Computing bamnado merged per-group-scaled heatmap matrix from bigWig files"

use rule heatmap_deeptools_per_sample_plot as heatmap_bamnado_merged_per_group_plot with:
    input:
        matrix=OUTPUT.select_heatmap_matrix(
            DataScalingTechnique.PER_GROUP,
            method=PileupMethod.BAMNADO,
            is_merged=True, scaling_method=_DEFAULT_SCALING_METHOD),
    output:
        heatmap=OUTPUT.select_heatmap_plot(
            DataScalingTechnique.PER_GROUP,
            method=PileupMethod.BAMNADO,
            is_merged=True, scaling_method=_DEFAULT_SCALING_METHOD),
    log: OUTPUT_DIR + "/logs/heatmap/merged/bamnado/" + DataScalingTechnique.PER_GROUP.value + "/heatmap.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/merged/bamnado/" + DataScalingTechnique.PER_GROUP.value + "/heatmap.tsv",
    message: "Generating bamnado merged per-group-scaled heatmap from matrix"

use rule heatmap_deeptools_per_sample_metaplot as heatmap_bamnado_merged_per_group_metaplot with:
    input:
        matrix=OUTPUT.select_heatmap_matrix(
            DataScalingTechnique.PER_GROUP,
            method=PileupMethod.BAMNADO,
            is_merged=True, scaling_method=_DEFAULT_SCALING_METHOD),
    output:
        metaplot=OUTPUT.select_heatmap_metaplot(
            DataScalingTechnique.PER_GROUP,
            method=PileupMethod.BAMNADO,
            is_merged=True, scaling_method=_DEFAULT_SCALING_METHOD),
    log: OUTPUT_DIR + "/logs/heatmap/merged/bamnado/" + DataScalingTechnique.PER_GROUP.value + "/metaplot.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/merged/bamnado/" + DataScalingTechnique.PER_GROUP.value + "/metaplot.tsv",
    message: "Generating bamnado merged per-group-scaled metaplot from matrix"


# ─── bamnado · merged · spikein · orlando ────────────────────────────────────
if _has_spikein_orlando:
    use rule heatmap_deeptools_per_sample_matrix as heatmap_bamnado_merged_spikein_orlando_matrix with:
        input:
            bigwigs=OUTPUT.select_bigwig_subtype(
                method=PileupMethod.BAMNADO,
                scale=DataScalingTechnique.SPIKEIN,
                spikein_method=SpikeInMethod.ORLANDO.value,
                is_merged=True,
            ),
        output:
            matrix=temp(OUTPUT.select_heatmap_matrix(
            DataScalingTechnique.SPIKEIN, 
            method=PileupMethod.BAMNADO,
            is_merged=True,
            spikein_method=SpikeInMethod.ORLANDO.value)),
        log: OUTPUT_DIR + "/logs/heatmap/merged/bamnado/spikein/orlando/matrix.log",
        benchmark: OUTPUT_DIR + "/.benchmark/heatmap/merged/bamnado/spikein/orlando/matrix.tsv",
        message: "Computing bamnado merged orlando spike-in heatmap matrix from bigWig files"

    use rule heatmap_deeptools_per_sample_plot as heatmap_bamnado_merged_spikein_orlando_plot with:
        input:
            matrix=OUTPUT.select_heatmap_matrix(
            DataScalingTechnique.SPIKEIN, 
            method=PileupMethod.BAMNADO,
            is_merged=True,
            spikein_method=SpikeInMethod.ORLANDO.value),
        output:
            heatmap=OUTPUT.select_heatmap_plot(
            DataScalingTechnique.SPIKEIN, 
            method=PileupMethod.BAMNADO,
            is_merged=True,
            spikein_method=SpikeInMethod.ORLANDO.value),
        log: OUTPUT_DIR + "/logs/heatmap/merged/bamnado/spikein/orlando/heatmap.log",
        benchmark: OUTPUT_DIR + "/.benchmark/heatmap/merged/bamnado/spikein/orlando/heatmap.tsv",
        message: "Generating bamnado merged orlando spike-in heatmap from matrix"

    use rule heatmap_deeptools_per_sample_metaplot as heatmap_bamnado_merged_spikein_orlando_metaplot with:
        input:
            matrix=OUTPUT.select_heatmap_matrix(
            DataScalingTechnique.SPIKEIN, 
            method=PileupMethod.BAMNADO,
            is_merged=True,
            spikein_method=SpikeInMethod.ORLANDO.value),
        output:
            metaplot=OUTPUT.select_heatmap_metaplot(
            DataScalingTechnique.SPIKEIN, 
            method=PileupMethod.BAMNADO,
            is_merged=True,
            spikein_method=SpikeInMethod.ORLANDO.value),
        log: OUTPUT_DIR + "/logs/heatmap/merged/bamnado/spikein/orlando/metaplot.log",
        benchmark: OUTPUT_DIR + "/.benchmark/heatmap/merged/bamnado/spikein/orlando/metaplot.tsv",
        message: "Generating bamnado merged orlando spike-in metaplot from matrix"


# ─── bamnado · merged · spikein · with_input ─────────────────────────────────
if _has_spikein_withinput:
    use rule heatmap_deeptools_per_sample_matrix as heatmap_bamnado_merged_spikein_withinput_matrix with:
        input:
            bigwigs=OUTPUT.select_bigwig_subtype(
                method=PileupMethod.BAMNADO,
                scale=DataScalingTechnique.SPIKEIN,
                spikein_method=SpikeInMethod.WITH_INPUT.value,
                is_merged=True,
            ),
        output:
            matrix=temp(OUTPUT.select_heatmap_matrix(
            DataScalingTechnique.SPIKEIN, 
            method=PileupMethod.BAMNADO,
            is_merged=True,
            spikein_method=SpikeInMethod.WITH_INPUT.value)),
        log: OUTPUT_DIR + "/logs/heatmap/merged/bamnado/spikein/with_input/matrix.log",
        benchmark: OUTPUT_DIR + "/.benchmark/heatmap/merged/bamnado/spikein/with_input/matrix.tsv",
        message: "Computing bamnado merged with_input spike-in heatmap matrix from bigWig files"

    use rule heatmap_deeptools_per_sample_plot as heatmap_bamnado_merged_spikein_withinput_plot with:
        input:
            matrix=OUTPUT.select_heatmap_matrix(
                DataScalingTechnique.SPIKEIN, 
                method=PileupMethod.BAMNADO, 
                is_merged=True, 
                spikein_method=SpikeInMethod.WITH_INPUT.value
            ),
        output:
            heatmap=OUTPUT.select_heatmap_plot(
            DataScalingTechnique.SPIKEIN, 
            method=PileupMethod.BAMNADO,
            is_merged=True,
            spikein_method=SpikeInMethod.WITH_INPUT.value),
        log: OUTPUT_DIR + "/logs/heatmap/merged/bamnado/spikein/with_input/heatmap.log",
        benchmark: OUTPUT_DIR + "/.benchmark/heatmap/merged/bamnado/spikein/with_input/heatmap.tsv",
        message: "Generating bamnado merged with_input spike-in heatmap from matrix"

    use rule heatmap_deeptools_per_sample_metaplot as heatmap_bamnado_merged_spikein_withinput_metaplot with:
        input:
            matrix=OUTPUT.select_heatmap_matrix(
            DataScalingTechnique.SPIKEIN, 
            method=PileupMethod.BAMNADO,
            is_merged=True,
            spikein_method=SpikeInMethod.WITH_INPUT.value),
        output:
            metaplot=OUTPUT.select_heatmap_metaplot(
            DataScalingTechnique.SPIKEIN, 
            method=PileupMethod.BAMNADO,
            is_merged=True,
            spikein_method=SpikeInMethod.WITH_INPUT.value),
        log: OUTPUT_DIR + "/logs/heatmap/merged/bamnado/spikein/with_input/metaplot.log",
        benchmark: OUTPUT_DIR + "/.benchmark/heatmap/merged/bamnado/spikein/with_input/metaplot.tsv",
        message: "Generating bamnado merged with_input spike-in metaplot from matrix"


# ─── homer · individual · per_sample ───────────────────────────────────────────
use rule heatmap_deeptools_per_sample_matrix as heatmap_homer_per_sample_matrix with:
    input:
        bigwigs=OUTPUT.select_bigwig_subtype(
            method=PileupMethod.HOMER,
            scale=DataScalingTechnique.PER_SAMPLE,
            ip_only=True,
        ),
    output:
        matrix=temp(OUTPUT.select_heatmap_matrix(
            DataScalingTechnique.PER_SAMPLE, 
            method=PileupMethod.HOMER)),
    log: OUTPUT_DIR + "/logs/heatmap/homer/" + DataScalingTechnique.PER_SAMPLE.value + "/matrix.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/homer/" + DataScalingTechnique.PER_SAMPLE.value + "/matrix.tsv",
    message: "Computing homer per_sample heatmap matrix from bigWig files"

use rule heatmap_deeptools_per_sample_plot as heatmap_homer_per_sample_plot with:
    input:
        matrix=OUTPUT.select_heatmap_matrix(
            DataScalingTechnique.PER_SAMPLE, 
            method=PileupMethod.HOMER),
    output:
        heatmap=OUTPUT.select_heatmap_plot(
            DataScalingTechnique.PER_SAMPLE, 
            method=PileupMethod.HOMER),
    log: OUTPUT_DIR + "/logs/heatmap/homer/" + DataScalingTechnique.PER_SAMPLE.value + "/heatmap.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/homer/" + DataScalingTechnique.PER_SAMPLE.value + "/heatmap.tsv",
    message: "Generating homer per_sample heatmap from matrix"

use rule heatmap_deeptools_per_sample_metaplot as heatmap_homer_per_sample_metaplot with:
    input:
        matrix=OUTPUT.select_heatmap_matrix(
            DataScalingTechnique.PER_SAMPLE, 
            method=PileupMethod.HOMER),
    output:
        metaplot=OUTPUT.select_heatmap_metaplot(
            DataScalingTechnique.PER_SAMPLE, 
            method=PileupMethod.HOMER),
    log: OUTPUT_DIR + "/logs/heatmap/homer/" + DataScalingTechnique.PER_SAMPLE.value + "/metaplot.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/homer/" + DataScalingTechnique.PER_SAMPLE.value + "/metaplot.tsv",
    message: "Generating homer per_sample metaplot from matrix"


# ─── homer · merged · per_sample ───────────────────────────────────────────────
use rule heatmap_deeptools_per_sample_matrix as heatmap_homer_merged_per_sample_matrix with:
    input:
        bigwigs=OUTPUT.select_bigwig_subtype(
            method=PileupMethod.HOMER,
            scale=DataScalingTechnique.PER_SAMPLE,
            is_merged=True,
        ),
    output:
        matrix=temp(OUTPUT.select_heatmap_matrix(
            DataScalingTechnique.PER_SAMPLE, 
            method=PileupMethod.HOMER,
            is_merged=True)),
    log: OUTPUT_DIR + "/logs/heatmap/merged/homer/" + DataScalingTechnique.PER_SAMPLE.value + "/matrix.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/merged/homer/" + DataScalingTechnique.PER_SAMPLE.value + "/matrix.tsv",
    message: "Computing homer merged per_sample heatmap matrix from bigWig files"

use rule heatmap_deeptools_per_sample_plot as heatmap_homer_merged_per_sample_plot with:
    input:
        matrix=OUTPUT.select_heatmap_matrix(
            DataScalingTechnique.PER_SAMPLE, 
            method=PileupMethod.HOMER,
            is_merged=True),
    output:
        heatmap=OUTPUT.select_heatmap_plot(
            DataScalingTechnique.PER_SAMPLE, 
            method=PileupMethod.HOMER,
            is_merged=True),
    log: OUTPUT_DIR + "/logs/heatmap/merged/homer/" + DataScalingTechnique.PER_SAMPLE.value + "/heatmap.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/merged/homer/" + DataScalingTechnique.PER_SAMPLE.value + "/heatmap.tsv",
    message: "Generating homer merged per_sample heatmap from matrix"

use rule heatmap_deeptools_per_sample_metaplot as heatmap_homer_merged_per_sample_metaplot with:
    input:
        matrix=OUTPUT.select_heatmap_matrix(
            DataScalingTechnique.PER_SAMPLE, 
            method=PileupMethod.HOMER, 
            is_merged=True
        ),
    output:
        metaplot=OUTPUT.select_heatmap_metaplot(
            DataScalingTechnique.PER_SAMPLE, 
            method=PileupMethod.HOMER, 
            is_merged=True
        ),
    log: OUTPUT_DIR + "/logs/heatmap/merged/homer/" + DataScalingTechnique.PER_SAMPLE.value + "/metaplot.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/merged/homer/" + DataScalingTechnique.PER_SAMPLE.value + "/metaplot.tsv",
    message: "Generating homer merged per_sample metaplot from matrix"
