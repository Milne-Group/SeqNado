from seqnado.workflow.helpers.common import define_memory_requested, define_time_requested
from seqnado import PileupMethod, DataScalingTechnique, SpikeInMethod

_spikein_cfg = getattr(CONFIG.assay_config, "spikein", None)
_spikein_methods = _spikein_cfg.method if _spikein_cfg else []
_has_spikein_orlando = SpikeInMethod.ORLANDO in _spikein_methods
_has_spikein_withinput = SpikeInMethod.WITH_INPUT in _spikein_methods

# Heatmap and metaplot generation from bigWig files.
# One set of rules per (pileup method) × (scale) × (merged) combination.
# Compatible combinations:
#   deeptools  : individual → unscaled/csaw/spikein_orlando/spikein_withinput ; merged → unscaled/csaw/spikein_orlando/spikein_withinput
#   bamnado    : individual → unscaled/csaw/spikein_orlando/spikein_withinput ; merged → unscaled/csaw/spikein_orlando/spikein_withinput
#   homer      : individual → unscaled              ; merged → unscaled


# ─── deeptools · individual · unscaled (base rules) ──────────────────────────
rule heatmap_deeptools_unscaled_matrix:
    input:
        bigwigs=OUTPUT.select_bigwig_subtype(
            method=PileupMethod.DEEPTOOLS,
            scale=DataScalingTechnique.UNSCALED,
            ip_only=True,
        ),
    output:
        matrix=temp(OUTPUT.select_heatmap_matrix(
            DataScalingTechnique.UNSCALED
            )),
    params:
        gtf=CONFIG.genome.gtf,
        options=str(CONFIG.third_party_tools.deeptools.compute_matrix.command_line_arguments),
    threads: CONFIG.third_party_tools.deeptools.compute_matrix.threads
    resources:
        runtime=lambda wildcards, attempt: f"{1 * 2**attempt}h",
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/heatmap/deeptools/unscaled/matrix.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/deeptools/unscaled/matrix.tsv",
    message: "Computing deeptools unscaled heatmap matrix from bigWig files"
    shell: """
    computeMatrix scale-regions -p {threads} {params.options} --smartLabels --missingDataAsZero -S {input.bigwigs} -R {params.gtf} -o {output.matrix} >> {log} 2>&1
    """


rule heatmap_deeptools_unscaled_plot:
    input:
        matrix=rules.heatmap_deeptools_unscaled_matrix.output.matrix,
    output:
        heatmap=OUTPUT.select_heatmap_plot(
            DataScalingTechnique.UNSCALED),
    params:
        options=str(CONFIG.third_party_tools.deeptools.plot_heatmap.command_line_arguments),
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/heatmap/deeptools/unscaled/heatmap.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/deeptools/unscaled/heatmap.tsv",
    message: "Generating deeptools unscaled heatmap from matrix"
    shell: """
    plotHeatmap -m {input.matrix} -out {output.heatmap} {params.options}
    """


rule heatmap_deeptools_unscaled_metaplot:
    input:
        matrix=rules.heatmap_deeptools_unscaled_matrix.output.matrix,
    output:
        metaplot=OUTPUT.select_heatmap_metaplot(
            DataScalingTechnique.UNSCALED),
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/heatmap/deeptools/unscaled/metaplot.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/deeptools/unscaled/metaplot.tsv",
    message: "Generating deeptools unscaled metaplot from heatmap matrix"
    shell: """
    plotProfile -m {input.matrix} -out {output.metaplot} --perGroup
    """


# ─── deeptools · individual · csaw ───────────────────────────────────────────
use rule heatmap_deeptools_unscaled_matrix as heatmap_deeptools_csaw_matrix with:
    input:
        bigwigs=OUTPUT.select_bigwig_subtype(
            method=PileupMethod.DEEPTOOLS,
            scale=DataScalingTechnique.CSAW,
            ip_only=True,
        ),
    output:
        matrix=temp(OUTPUT.select_heatmap_matrix(
            DataScalingTechnique.CSAW)),
    log: OUTPUT_DIR + "/logs/heatmap/deeptools/csaw/matrix.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/deeptools/csaw/matrix.tsv",
    message: "Computing deeptools CSAW-scaled heatmap matrix from bigWig files"

use rule heatmap_deeptools_unscaled_plot as heatmap_deeptools_csaw_plot with:
    input:
        matrix=OUTPUT.select_heatmap_matrix(
            DataScalingTechnique.CSAW),
    output:
        heatmap=OUTPUT.select_heatmap_plot(
            DataScalingTechnique.CSAW),
    log: OUTPUT_DIR + "/logs/heatmap/deeptools/csaw/heatmap.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/deeptools/csaw/heatmap.tsv",
    message: "Generating deeptools CSAW-scaled heatmap from matrix"

use rule heatmap_deeptools_unscaled_metaplot as heatmap_deeptools_csaw_metaplot with:
    input:
        matrix=OUTPUT.select_heatmap_matrix(
            DataScalingTechnique.CSAW),
    output:
        metaplot=OUTPUT.select_heatmap_metaplot(
            DataScalingTechnique.CSAW),
    log: OUTPUT_DIR + "/logs/heatmap/deeptools/csaw/metaplot.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/deeptools/csaw/metaplot.tsv",
    message: "Generating deeptools CSAW-scaled metaplot from matrix"


# ─── deeptools · individual · spikein · orlando ───────────────────────────────
if _has_spikein_orlando:
    use rule heatmap_deeptools_unscaled_matrix as heatmap_deeptools_spikein_orlando_matrix with:
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

    use rule heatmap_deeptools_unscaled_plot as heatmap_deeptools_spikein_orlando_plot with:
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

    use rule heatmap_deeptools_unscaled_metaplot as heatmap_deeptools_spikein_orlando_metaplot with:
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
    use rule heatmap_deeptools_unscaled_matrix as heatmap_deeptools_spikein_withinput_matrix with:
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

    use rule heatmap_deeptools_unscaled_plot as heatmap_deeptools_spikein_withinput_plot with:
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

    use rule heatmap_deeptools_unscaled_metaplot as heatmap_deeptools_spikein_withinput_metaplot with:
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


# ─── deeptools · merged · unscaled ───────────────────────────────────────────
use rule heatmap_deeptools_unscaled_matrix as heatmap_deeptools_merged_matrix with:
    input:
        bigwigs=OUTPUT.select_bigwig_subtype(
            method=PileupMethod.DEEPTOOLS,
            scale=DataScalingTechnique.UNSCALED,
            is_merged=True,
        ),
    output:
        matrix=temp(OUTPUT.select_heatmap_matrix(
            DataScalingTechnique.UNSCALED,
            is_merged=True)),
    log: OUTPUT_DIR + "/logs/heatmap/merged/deeptools/unscaled/matrix.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/merged/deeptools/unscaled/matrix.tsv",
    message: "Computing deeptools merged unscaled heatmap matrix from bigWig files"

use rule heatmap_deeptools_unscaled_plot as heatmap_deeptools_merged_plot with:
    input:
        matrix=OUTPUT.select_heatmap_matrix(
            DataScalingTechnique.UNSCALED,
            is_merged=True),
    output:
        heatmap=OUTPUT.select_heatmap_plot(
            DataScalingTechnique.UNSCALED,
            is_merged=True),
    log: OUTPUT_DIR + "/logs/heatmap/merged/deeptools/unscaled/heatmap.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/merged/deeptools/unscaled/heatmap.tsv",
    message: "Generating deeptools merged unscaled heatmap from matrix"

use rule heatmap_deeptools_unscaled_metaplot as heatmap_deeptools_merged_metaplot with:
    input:
        matrix=OUTPUT.select_heatmap_matrix(
            DataScalingTechnique.UNSCALED,
            is_merged=True),
    output:
        metaplot=OUTPUT.select_heatmap_metaplot(
            DataScalingTechnique.UNSCALED,
            is_merged=True),
    log: OUTPUT_DIR + "/logs/heatmap/merged/deeptools/unscaled/metaplot.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/merged/deeptools/unscaled/metaplot.tsv",
    message: "Generating deeptools merged unscaled metaplot from matrix"


# ─── deeptools · merged · csaw ───────────────────────────────────────────────
use rule heatmap_deeptools_unscaled_matrix as heatmap_deeptools_merged_csaw_matrix with:
    input:
        bigwigs=OUTPUT.select_bigwig_subtype(
            method=PileupMethod.DEEPTOOLS,
            scale=DataScalingTechnique.CSAW,
            is_merged=True,
        ),
    output:
        matrix=temp(OUTPUT.select_heatmap_matrix(
            DataScalingTechnique.CSAW,
            is_merged=True)),
    log: OUTPUT_DIR + "/logs/heatmap/merged/deeptools/csaw/matrix.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/merged/deeptools/csaw/matrix.tsv",
    message: "Computing deeptools merged CSAW-scaled heatmap matrix from bigWig files"

use rule heatmap_deeptools_unscaled_plot as heatmap_deeptools_merged_csaw_plot with:
    input:
        matrix=OUTPUT.select_heatmap_matrix(
            DataScalingTechnique.CSAW,
            is_merged=True),
    output:
        heatmap=OUTPUT.select_heatmap_plot(
            DataScalingTechnique.CSAW,
            is_merged=True),
    log: OUTPUT_DIR + "/logs/heatmap/merged/deeptools/csaw/heatmap.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/merged/deeptools/csaw/heatmap.tsv",
    message: "Generating deeptools merged CSAW-scaled heatmap from matrix"

use rule heatmap_deeptools_unscaled_metaplot as heatmap_deeptools_merged_csaw_metaplot with:
    input:
        matrix=OUTPUT.select_heatmap_matrix(
            DataScalingTechnique.CSAW,
            is_merged=True),
    output:
        metaplot=OUTPUT.select_heatmap_metaplot(
            DataScalingTechnique.CSAW,
            is_merged=True),
    log: OUTPUT_DIR + "/logs/heatmap/merged/deeptools/csaw/metaplot.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/merged/deeptools/csaw/metaplot.tsv",
    message: "Generating deeptools merged CSAW-scaled metaplot from matrix"


# ─── deeptools · merged · spikein · orlando ──────────────────────────────────
if _has_spikein_orlando:
    use rule heatmap_deeptools_unscaled_matrix as heatmap_deeptools_merged_spikein_orlando_matrix with:
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

    use rule heatmap_deeptools_unscaled_plot as heatmap_deeptools_merged_spikein_orlando_plot with:
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

    use rule heatmap_deeptools_unscaled_metaplot as heatmap_deeptools_merged_spikein_orlando_metaplot with:
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
    use rule heatmap_deeptools_unscaled_matrix as heatmap_deeptools_merged_spikein_withinput_matrix with:
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

    use rule heatmap_deeptools_unscaled_plot as heatmap_deeptools_merged_spikein_withinput_plot with:
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

    use rule heatmap_deeptools_unscaled_metaplot as heatmap_deeptools_merged_spikein_withinput_metaplot with:
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


# ─── bamnado · individual · unscaled ─────────────────────────────────────────
use rule heatmap_deeptools_unscaled_matrix as heatmap_bamnado_unscaled_matrix with:
    input:
        bigwigs=OUTPUT.select_bigwig_subtype(
            method=PileupMethod.BAMNADO,
            scale=DataScalingTechnique.UNSCALED,
            ip_only=True,
        ),
    output:
        matrix=temp(OUTPUT.select_heatmap_matrix(
            DataScalingTechnique.UNSCALED, 
            method=PileupMethod.BAMNADO)),
    log: OUTPUT_DIR + "/logs/heatmap/bamnado/unscaled/matrix.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/bamnado/unscaled/matrix.tsv",
    message: "Computing bamnado unscaled heatmap matrix from bigWig files"

use rule heatmap_deeptools_unscaled_plot as heatmap_bamnado_unscaled_plot with:
    input:
        matrix=OUTPUT.select_heatmap_matrix(
            DataScalingTechnique.UNSCALED, 
            method=PileupMethod.BAMNADO),
    output:
        heatmap=OUTPUT.select_heatmap_plot(
            DataScalingTechnique.UNSCALED, 
            method=PileupMethod.BAMNADO),
    log: OUTPUT_DIR + "/logs/heatmap/bamnado/unscaled/heatmap.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/bamnado/unscaled/heatmap.tsv",
    message: "Generating bamnado unscaled heatmap from matrix"

use rule heatmap_deeptools_unscaled_metaplot as heatmap_bamnado_unscaled_metaplot with:
    input:
        matrix=OUTPUT.select_heatmap_matrix(
            DataScalingTechnique.UNSCALED, 
            method=PileupMethod.BAMNADO),
    output:
        metaplot=OUTPUT.select_heatmap_metaplot(
            DataScalingTechnique.UNSCALED, 
            method=PileupMethod.BAMNADO),
    log: OUTPUT_DIR + "/logs/heatmap/bamnado/unscaled/metaplot.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/bamnado/unscaled/metaplot.tsv",
    message: "Generating bamnado unscaled metaplot from matrix"


# ─── bamnado · merged · unscaled ─────────────────────────────────────────────
use rule heatmap_deeptools_unscaled_matrix as heatmap_bamnado_merged_unscaled_matrix with:
    input:
        bigwigs=OUTPUT.select_bigwig_subtype(
            method=PileupMethod.BAMNADO,
            scale=DataScalingTechnique.UNSCALED,
            is_merged=True,
        ),
    output:
        matrix=temp(OUTPUT.select_heatmap_matrix(
            DataScalingTechnique.UNSCALED, 
            method=PileupMethod.BAMNADO,
            is_merged=True)),
    log: OUTPUT_DIR + "/logs/heatmap/merged/bamnado/unscaled/matrix.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/merged/bamnado/unscaled/matrix.tsv",
    message: "Computing bamnado merged unscaled heatmap matrix from bigWig files"

use rule heatmap_deeptools_unscaled_plot as heatmap_bamnado_merged_unscaled_plot with:
    input:
        matrix=OUTPUT.select_heatmap_matrix(
            DataScalingTechnique.UNSCALED, 
            method=PileupMethod.BAMNADO,
            is_merged=True),
    output:
        heatmap=OUTPUT.select_heatmap_plot(
            DataScalingTechnique.UNSCALED, 
            method=PileupMethod.BAMNADO,
            is_merged=True),
    log: OUTPUT_DIR + "/logs/heatmap/merged/bamnado/unscaled/heatmap.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/merged/bamnado/unscaled/heatmap.tsv",
    message: "Generating bamnado merged unscaled heatmap from matrix"

use rule heatmap_deeptools_unscaled_metaplot as heatmap_bamnado_merged_unscaled_metaplot with:
    input:
        matrix=OUTPUT.select_heatmap_matrix(
            DataScalingTechnique.UNSCALED, 
            method=PileupMethod.BAMNADO,
            is_merged=True),
    output:
        metaplot=OUTPUT.select_heatmap_metaplot(
            DataScalingTechnique.UNSCALED, 
            method=PileupMethod.BAMNADO,
            is_merged=True),
    log: OUTPUT_DIR + "/logs/heatmap/merged/bamnado/unscaled/metaplot.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/merged/bamnado/unscaled/metaplot.tsv",
    message: "Generating bamnado merged unscaled metaplot from matrix"


# ─── bamnado · merged · csaw ─────────────────────────────────────────────────
use rule heatmap_deeptools_unscaled_matrix as heatmap_bamnado_merged_csaw_matrix with:
    input:
        bigwigs=OUTPUT.select_bigwig_subtype(
            method=PileupMethod.BAMNADO,
            scale=DataScalingTechnique.CSAW,
            is_merged=True,
        ),
    output:
        matrix=temp(OUTPUT.select_heatmap_matrix(
            DataScalingTechnique.CSAW, 
            method=PileupMethod.BAMNADO,
            is_merged=True)),
    log: OUTPUT_DIR + "/logs/heatmap/merged/bamnado/csaw/matrix.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/merged/bamnado/csaw/matrix.tsv",
    message: "Computing bamnado merged CSAW-scaled heatmap matrix from bigWig files"

use rule heatmap_deeptools_unscaled_plot as heatmap_bamnado_merged_csaw_plot with:
    input:
        matrix=OUTPUT.select_heatmap_matrix(
            DataScalingTechnique.CSAW, 
            method=PileupMethod.BAMNADO,
            is_merged=True),
    output:
        heatmap=OUTPUT.select_heatmap_plot(
            DataScalingTechnique.CSAW, 
            method=PileupMethod.BAMNADO,
            is_merged=True),
    log: OUTPUT_DIR + "/logs/heatmap/merged/bamnado/csaw/heatmap.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/merged/bamnado/csaw/heatmap.tsv",
    message: "Generating bamnado merged CSAW-scaled heatmap from matrix"

use rule heatmap_deeptools_unscaled_metaplot as heatmap_bamnado_merged_csaw_metaplot with:
    input:
        matrix=OUTPUT.select_heatmap_matrix(
            DataScalingTechnique.CSAW, 
            method=PileupMethod.BAMNADO,
            is_merged=True),
    output:
        metaplot=OUTPUT.select_heatmap_metaplot(
            DataScalingTechnique.CSAW, 
            method=PileupMethod.BAMNADO,
            is_merged=True),
    log: OUTPUT_DIR + "/logs/heatmap/merged/bamnado/csaw/metaplot.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/merged/bamnado/csaw/metaplot.tsv",
    message: "Generating bamnado merged CSAW-scaled metaplot from matrix"


# ─── bamnado · merged · spikein · orlando ────────────────────────────────────
if _has_spikein_orlando:
    use rule heatmap_deeptools_unscaled_matrix as heatmap_bamnado_merged_spikein_orlando_matrix with:
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

    use rule heatmap_deeptools_unscaled_plot as heatmap_bamnado_merged_spikein_orlando_plot with:
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

    use rule heatmap_deeptools_unscaled_metaplot as heatmap_bamnado_merged_spikein_orlando_metaplot with:
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
    use rule heatmap_deeptools_unscaled_matrix as heatmap_bamnado_merged_spikein_withinput_matrix with:
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

    use rule heatmap_deeptools_unscaled_plot as heatmap_bamnado_merged_spikein_withinput_plot with:
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

    use rule heatmap_deeptools_unscaled_metaplot as heatmap_bamnado_merged_spikein_withinput_metaplot with:
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


# ─── homer · individual · unscaled ───────────────────────────────────────────
use rule heatmap_deeptools_unscaled_matrix as heatmap_homer_unscaled_matrix with:
    input:
        bigwigs=OUTPUT.select_bigwig_subtype(
            method=PileupMethod.HOMER,
            scale=DataScalingTechnique.UNSCALED,
            ip_only=True,
        ),
    output:
        matrix=temp(OUTPUT.select_heatmap_matrix(
            DataScalingTechnique.UNSCALED, 
            method=PileupMethod.HOMER)),
    log: OUTPUT_DIR + "/logs/heatmap/homer/unscaled/matrix.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/homer/unscaled/matrix.tsv",
    message: "Computing homer unscaled heatmap matrix from bigWig files"

use rule heatmap_deeptools_unscaled_plot as heatmap_homer_unscaled_plot with:
    input:
        matrix=OUTPUT.select_heatmap_matrix(
            DataScalingTechnique.UNSCALED, 
            method=PileupMethod.HOMER),
    output:
        heatmap=OUTPUT.select_heatmap_plot(
            DataScalingTechnique.UNSCALED, 
            method=PileupMethod.HOMER),
    log: OUTPUT_DIR + "/logs/heatmap/homer/unscaled/heatmap.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/homer/unscaled/heatmap.tsv",
    message: "Generating homer unscaled heatmap from matrix"

use rule heatmap_deeptools_unscaled_metaplot as heatmap_homer_unscaled_metaplot with:
    input:
        matrix=OUTPUT.select_heatmap_matrix(
            DataScalingTechnique.UNSCALED, 
            method=PileupMethod.HOMER),
    output:
        metaplot=OUTPUT.select_heatmap_metaplot(
            DataScalingTechnique.UNSCALED, 
            method=PileupMethod.HOMER),
    log: OUTPUT_DIR + "/logs/heatmap/homer/unscaled/metaplot.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/homer/unscaled/metaplot.tsv",
    message: "Generating homer unscaled metaplot from matrix"


# ─── homer · merged · unscaled ───────────────────────────────────────────────
use rule heatmap_deeptools_unscaled_matrix as heatmap_homer_merged_unscaled_matrix with:
    input:
        bigwigs=OUTPUT.select_bigwig_subtype(
            method=PileupMethod.HOMER,
            scale=DataScalingTechnique.UNSCALED,
            is_merged=True,
        ),
    output:
        matrix=temp(OUTPUT.select_heatmap_matrix(
            DataScalingTechnique.UNSCALED, 
            method=PileupMethod.HOMER,
            is_merged=True)),
    log: OUTPUT_DIR + "/logs/heatmap/merged/homer/unscaled/matrix.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/merged/homer/unscaled/matrix.tsv",
    message: "Computing homer merged unscaled heatmap matrix from bigWig files"

use rule heatmap_deeptools_unscaled_plot as heatmap_homer_merged_unscaled_plot with:
    input:
        matrix=OUTPUT.select_heatmap_matrix(
            DataScalingTechnique.UNSCALED, 
            method=PileupMethod.HOMER,
            is_merged=True),
    output:
        heatmap=OUTPUT.select_heatmap_plot(
            DataScalingTechnique.UNSCALED, 
            method=PileupMethod.HOMER,
            is_merged=True),
    log: OUTPUT_DIR + "/logs/heatmap/merged/homer/unscaled/heatmap.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/merged/homer/unscaled/heatmap.tsv",
    message: "Generating homer merged unscaled heatmap from matrix"

use rule heatmap_deeptools_unscaled_metaplot as heatmap_homer_merged_unscaled_metaplot with:
    input:
        matrix=OUTPUT.select_heatmap_matrix(
            DataScalingTechnique.UNSCALED, 
            method=PileupMethod.HOMER, 
            is_merged=True
        ),
    output:
        metaplot=OUTPUT.select_heatmap_metaplot(
            DataScalingTechnique.UNSCALED, 
            method=PileupMethod.HOMER, 
            is_merged=True
        ),
    log: OUTPUT_DIR + "/logs/heatmap/merged/homer/unscaled/metaplot.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/merged/homer/unscaled/metaplot.tsv",
    message: "Generating homer merged unscaled metaplot from matrix"
