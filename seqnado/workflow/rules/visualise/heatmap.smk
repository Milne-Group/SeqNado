from seqnado.workflow.helpers.common import define_memory_requested, define_time_requested
from seqnado import PileupMethod, DataScalingTechnique

# Heatmap generation from bigWig files for scaled, unscaled, spike-in
# and merged unscaled, scaled and spike-in data.
# with both heatmap and metaplot outputs.
# Uses deeptools computeMatrix and plotHeatmap/plotProfile.

rule heatmap_matrix:
    input:
        bigwigs=OUTPUT.select_bigwig_subtype(
            method=PileupMethod.DEEPTOOLS,
            scale=DataScalingTechnique.UNSCALED
        ),
    output:
        matrix=temp(OUTPUT_DIR + "/heatmap/unscaled/matrix.mat.gz"),
    params:
        gtf=CONFIG.genome.gtf,
        options=str(CONFIG.third_party_tools.deeptools.compute_matrix.command_line_arguments),
    threads: CONFIG.third_party_tools.deeptools.compute_matrix.threads
    resources:
        runtime=lambda wildcards, attempt: f"{1 * 2**attempt}h",
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/heatmap/unscaled/matrix.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/unscaled/matrix.tsv",
    message: "Computing heatmap matrix from bigWig files"
    shell: """
    computeMatrix scale-regions -p {threads} {params.options} --smartLabels --missingDataAsZero -S {input.bigwigs} -R {params.gtf} -o {output.matrix} >> {log} 2>&1
    """


rule heatmap_plot:
    input:
        matrix=OUTPUT_DIR + "/heatmap/unscaled/matrix.mat.gz",
    output:
        heatmap=OUTPUT_DIR + "/heatmap/unscaled/heatmap.pdf",
    params:
        options=str(CONFIG.third_party_tools.deeptools.plot_heatmap.command_line_arguments),
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/heatmap/unscaled/heatmap.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/unscaled/heatmap.tsv",
    message: "Generating heatmap from matrix"
    shell: """
    plotHeatmap -m {input.matrix} -out {output.heatmap} {params.options}
    """

rule heatmap_metaplot:
    input:
        matrix=OUTPUT_DIR + "/heatmap/unscaled/matrix.mat.gz",
    output:
        metaplot=OUTPUT_DIR + "/heatmap/unscaled/metaplot.pdf",
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/heatmap/unscaled/metaplot.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/unscaled/metaplot.tsv",
    message: "Generating metaplot from heatmap matrix"
    shell: """
    plotProfile -m {input.matrix} -out {output.metaplot} --perGroup
    """


# --- CSAW scaled ---

use rule heatmap_matrix as heatmap_matrix_csaw with:
    input:
        bigwigs=OUTPUT.select_bigwig_subtype(
            method=PileupMethod.DEEPTOOLS,
            scale=DataScalingTechnique.CSAW
        ),
    output:
        matrix=temp(OUTPUT_DIR + "/heatmap/csaw/matrix.mat.gz"),
    log: OUTPUT_DIR + "/logs/heatmap/csaw/matrix.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/csaw/matrix.tsv",
    message: "Computing CSAW-scaled heatmap matrix from bigWig files"

use rule heatmap_plot as heatmap_plot_csaw with:
    input:
        matrix=OUTPUT_DIR + "/heatmap/csaw/matrix.mat.gz",
    output:
        heatmap=OUTPUT_DIR + "/heatmap/csaw/heatmap.pdf",
    log: OUTPUT_DIR + "/logs/heatmap/csaw/heatmap.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/csaw/heatmap.tsv",
    message: "Generating CSAW-scaled heatmap from matrix"

use rule heatmap_metaplot as heatmap_metaplot_csaw with:
    input:
        matrix=OUTPUT_DIR + "/heatmap/csaw/matrix.mat.gz",
    output:
        metaplot=OUTPUT_DIR + "/heatmap/csaw/metaplot.pdf",
    log: OUTPUT_DIR + "/logs/heatmap/csaw/metaplot.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/csaw/metaplot.tsv",
    message: "Generating CSAW-scaled metaplot from matrix"


# --- Spike-in normalized ---

use rule heatmap_matrix as heatmap_matrix_spikein with:
    input:
        bigwigs=OUTPUT.select_bigwig_subtype(
            method=PileupMethod.DEEPTOOLS,
            scale=DataScalingTechnique.SPIKEIN
        ),
    output:
        matrix=temp(OUTPUT_DIR + "/heatmap/spikein/matrix.mat.gz"),
    log: OUTPUT_DIR + "/logs/heatmap/spikein/matrix.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/spikein/matrix.tsv",
    message: "Computing spike-in normalized heatmap matrix from bigWig files"

use rule heatmap_plot as heatmap_plot_spikein with:
    input:
        matrix=OUTPUT_DIR + "/heatmap/spikein/matrix.mat.gz",
    output:
        heatmap=OUTPUT_DIR + "/heatmap/spikein/heatmap.pdf",
    log: OUTPUT_DIR + "/logs/heatmap/spikein/heatmap.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/spikein/heatmap.tsv",
    message: "Generating spike-in normalized heatmap from matrix"

use rule heatmap_metaplot as heatmap_metaplot_spikein with:
    input:
        matrix=OUTPUT_DIR + "/heatmap/spikein/matrix.mat.gz",
    output:
        metaplot=OUTPUT_DIR + "/heatmap/spikein/metaplot.pdf",
    log: OUTPUT_DIR + "/logs/heatmap/spikein/metaplot.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/spikein/metaplot.tsv",
    message: "Generating spike-in normalized metaplot from matrix"


