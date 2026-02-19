from seqnado.workflow.helpers.common import define_memory_requested, define_time_requested
from seqnado import PileupMethod, DataScalingTechnique

# Heatmap and metaplot generation from bigWig files.
# One set of rules per (pileup method) × (scale) × (merged) combination.
# Compatible combinations:
#   deeptools  : individual → unscaled/csaw/spikein ; merged → unscaled/csaw/spikein
#   bamnado    : individual → unscaled              ; merged → unscaled/csaw/spikein
#   homer      : individual → unscaled              ; merged → unscaled


# ─── deeptools · individual · unscaled (base rules) ──────────────────────────

rule heatmap_matrix:
    input:
        bigwigs=OUTPUT.select_bigwig_subtype(
            method=PileupMethod.DEEPTOOLS,
            scale=DataScalingTechnique.UNSCALED,
            ip_only=True,
        ),
    output:
        matrix=temp(OUTPUT.select_heatmap_matrix(DataScalingTechnique.UNSCALED)),
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


rule heatmap_plot:
    input:
        matrix=rules.heatmap_matrix.output.matrix,
    output:
        heatmap=OUTPUT.select_heatmap_plot(DataScalingTechnique.UNSCALED),
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


rule heatmap_metaplot:
    input:
        matrix=OUTPUT.select_heatmap_matrix(DataScalingTechnique.UNSCALED),
    output:
        metaplot=OUTPUT.select_heatmap_metaplot(DataScalingTechnique.UNSCALED),
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

use rule heatmap_matrix as heatmap_matrix_csaw with:
    input:
        bigwigs=OUTPUT.select_bigwig_subtype(
            method=PileupMethod.DEEPTOOLS,
            scale=DataScalingTechnique.CSAW,
            ip_only=True,
        ),
    output:
        matrix=temp(OUTPUT.select_heatmap_matrix(DataScalingTechnique.CSAW)),
    log: OUTPUT_DIR + "/logs/heatmap/deeptools/csaw/matrix.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/deeptools/csaw/matrix.tsv",
    message: "Computing deeptools CSAW-scaled heatmap matrix from bigWig files"

use rule heatmap_plot as heatmap_plot_csaw with:
    input:
        matrix=OUTPUT.select_heatmap_matrix(DataScalingTechnique.CSAW),
    output:
        heatmap=OUTPUT.select_heatmap_plot(DataScalingTechnique.CSAW),
    log: OUTPUT_DIR + "/logs/heatmap/deeptools/csaw/heatmap.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/deeptools/csaw/heatmap.tsv",
    message: "Generating deeptools CSAW-scaled heatmap from matrix"

use rule heatmap_metaplot as heatmap_metaplot_csaw with:
    input:
        matrix=OUTPUT.select_heatmap_matrix(DataScalingTechnique.CSAW),
    output:
        metaplot=OUTPUT.select_heatmap_metaplot(DataScalingTechnique.CSAW),
    log: OUTPUT_DIR + "/logs/heatmap/deeptools/csaw/metaplot.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/deeptools/csaw/metaplot.tsv",
    message: "Generating deeptools CSAW-scaled metaplot from matrix"


# ─── deeptools · individual · spikein ────────────────────────────────────────

use rule heatmap_matrix as heatmap_matrix_spikein with:
    input:
        bigwigs=OUTPUT.select_bigwig_subtype(
            method=PileupMethod.DEEPTOOLS,
            scale=DataScalingTechnique.SPIKEIN,
            ip_only=True,
        ),
    output:
        matrix=temp(OUTPUT.select_heatmap_matrix(DataScalingTechnique.SPIKEIN)),
    log: OUTPUT_DIR + "/logs/heatmap/deeptools/spikein/matrix.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/deeptools/spikein/matrix.tsv",
    message: "Computing deeptools spike-in normalized heatmap matrix from bigWig files"

use rule heatmap_plot as heatmap_plot_spikein with:
    input:
        matrix=OUTPUT.select_heatmap_matrix(DataScalingTechnique.SPIKEIN),
    output:
        heatmap=OUTPUT.select_heatmap_plot(DataScalingTechnique.SPIKEIN),
    log: OUTPUT_DIR + "/logs/heatmap/deeptools/spikein/heatmap.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/deeptools/spikein/heatmap.tsv",
    message: "Generating deeptools spike-in normalized heatmap from matrix"

use rule heatmap_metaplot as heatmap_metaplot_spikein with:
    input:
        matrix=OUTPUT.select_heatmap_matrix(DataScalingTechnique.SPIKEIN),
    output:
        metaplot=OUTPUT.select_heatmap_metaplot(DataScalingTechnique.SPIKEIN),
    log: OUTPUT_DIR + "/logs/heatmap/deeptools/spikein/metaplot.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/deeptools/spikein/metaplot.tsv",
    message: "Generating deeptools spike-in normalized metaplot from matrix"


# ─── deeptools · merged · unscaled ───────────────────────────────────────────

use rule heatmap_matrix as heatmap_matrix_merged with:
    input:
        bigwigs=OUTPUT.select_bigwig_subtype(
            method=PileupMethod.DEEPTOOLS,
            scale=DataScalingTechnique.UNSCALED,
            is_merged=True,
        ),
    output:
        matrix=temp(OUTPUT.select_heatmap_matrix(DataScalingTechnique.UNSCALED, is_merged=True)),
    log: OUTPUT_DIR + "/logs/heatmap/merged/deeptools/unscaled/matrix.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/merged/deeptools/unscaled/matrix.tsv",
    message: "Computing deeptools merged unscaled heatmap matrix from bigWig files"

use rule heatmap_plot as heatmap_plot_merged with:
    input:
        matrix=OUTPUT.select_heatmap_matrix(DataScalingTechnique.UNSCALED, is_merged=True),
    output:
        heatmap=OUTPUT.select_heatmap_plot(DataScalingTechnique.UNSCALED, is_merged=True),
    log: OUTPUT_DIR + "/logs/heatmap/merged/deeptools/unscaled/heatmap.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/merged/deeptools/unscaled/heatmap.tsv",
    message: "Generating deeptools merged unscaled heatmap from matrix"

use rule heatmap_metaplot as heatmap_metaplot_merged with:
    input:
        matrix=OUTPUT.select_heatmap_matrix(DataScalingTechnique.UNSCALED, is_merged=True),
    output:
        metaplot=OUTPUT.select_heatmap_metaplot(DataScalingTechnique.UNSCALED, is_merged=True),
    log: OUTPUT_DIR + "/logs/heatmap/merged/deeptools/unscaled/metaplot.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/merged/deeptools/unscaled/metaplot.tsv",
    message: "Generating deeptools merged unscaled metaplot from matrix"


# ─── deeptools · merged · csaw ───────────────────────────────────────────────

use rule heatmap_matrix as heatmap_matrix_merged_csaw with:
    input:
        bigwigs=OUTPUT.select_bigwig_subtype(
            method=PileupMethod.DEEPTOOLS,
            scale=DataScalingTechnique.CSAW,
            is_merged=True,
        ),
    output:
        matrix=temp(OUTPUT.select_heatmap_matrix(DataScalingTechnique.CSAW, is_merged=True)),
    log: OUTPUT_DIR + "/logs/heatmap/merged/deeptools/csaw/matrix.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/merged/deeptools/csaw/matrix.tsv",
    message: "Computing deeptools merged CSAW-scaled heatmap matrix from bigWig files"

use rule heatmap_plot as heatmap_plot_merged_csaw with:
    input:
        matrix=OUTPUT.select_heatmap_matrix(DataScalingTechnique.CSAW, is_merged=True),
    output:
        heatmap=OUTPUT.select_heatmap_plot(DataScalingTechnique.CSAW, is_merged=True),
    log: OUTPUT_DIR + "/logs/heatmap/merged/deeptools/csaw/heatmap.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/merged/deeptools/csaw/heatmap.tsv",
    message: "Generating deeptools merged CSAW-scaled heatmap from matrix"

use rule heatmap_metaplot as heatmap_metaplot_merged_csaw with:
    input:
        matrix=OUTPUT.select_heatmap_matrix(DataScalingTechnique.CSAW, is_merged=True),
    output:
        metaplot=OUTPUT.select_heatmap_metaplot(DataScalingTechnique.CSAW, is_merged=True),
    log: OUTPUT_DIR + "/logs/heatmap/merged/deeptools/csaw/metaplot.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/merged/deeptools/csaw/metaplot.tsv",
    message: "Generating deeptools merged CSAW-scaled metaplot from matrix"


# ─── deeptools · merged · spikein ────────────────────────────────────────────

use rule heatmap_matrix as heatmap_matrix_merged_spikein with:
    input:
        bigwigs=OUTPUT.select_bigwig_subtype(
            method=PileupMethod.DEEPTOOLS,
            scale=DataScalingTechnique.SPIKEIN,
            is_merged=True,
        ),
    output:
        matrix=temp(OUTPUT.select_heatmap_matrix(DataScalingTechnique.SPIKEIN, is_merged=True)),
    log: OUTPUT_DIR + "/logs/heatmap/merged/deeptools/spikein/matrix.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/merged/deeptools/spikein/matrix.tsv",
    message: "Computing deeptools merged spike-in normalized heatmap matrix from bigWig files"

use rule heatmap_plot as heatmap_plot_merged_spikein with:
    input:
        matrix=OUTPUT.select_heatmap_matrix(DataScalingTechnique.SPIKEIN, is_merged=True),
    output:
        heatmap=OUTPUT.select_heatmap_plot(DataScalingTechnique.SPIKEIN, is_merged=True),
    log: OUTPUT_DIR + "/logs/heatmap/merged/deeptools/spikein/heatmap.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/merged/deeptools/spikein/heatmap.tsv",
    message: "Generating deeptools merged spike-in normalized heatmap from matrix"

use rule heatmap_metaplot as heatmap_metaplot_merged_spikein with:
    input:
        matrix=OUTPUT.select_heatmap_matrix(DataScalingTechnique.SPIKEIN, is_merged=True),
    output:
        metaplot=OUTPUT.select_heatmap_metaplot(DataScalingTechnique.SPIKEIN, is_merged=True),
    log: OUTPUT_DIR + "/logs/heatmap/merged/deeptools/spikein/metaplot.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/merged/deeptools/spikein/metaplot.tsv",
    message: "Generating deeptools merged spike-in normalized metaplot from matrix"


# ─── bamnado · individual · unscaled ─────────────────────────────────────────

use rule heatmap_matrix as heatmap_matrix_bamnado with:
    input:
        bigwigs=OUTPUT.select_bigwig_subtype(
            method=PileupMethod.BAMNADO,
            scale=DataScalingTechnique.UNSCALED,
            ip_only=True,
        ),
    output:
        matrix=temp(OUTPUT.select_heatmap_matrix(DataScalingTechnique.UNSCALED, method=PileupMethod.BAMNADO)),
    log: OUTPUT_DIR + "/logs/heatmap/bamnado/unscaled/matrix.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/bamnado/unscaled/matrix.tsv",
    message: "Computing bamnado unscaled heatmap matrix from bigWig files"

use rule heatmap_plot as heatmap_plot_bamnado with:
    input:
        matrix=OUTPUT.select_heatmap_matrix(DataScalingTechnique.UNSCALED, method=PileupMethod.BAMNADO),
    output:
        heatmap=OUTPUT.select_heatmap_plot(DataScalingTechnique.UNSCALED, method=PileupMethod.BAMNADO),
    log: OUTPUT_DIR + "/logs/heatmap/bamnado/unscaled/heatmap.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/bamnado/unscaled/heatmap.tsv",
    message: "Generating bamnado unscaled heatmap from matrix"

use rule heatmap_metaplot as heatmap_metaplot_bamnado with:
    input:
        matrix=OUTPUT.select_heatmap_matrix(DataScalingTechnique.UNSCALED, method=PileupMethod.BAMNADO),
    output:
        metaplot=OUTPUT.select_heatmap_metaplot(DataScalingTechnique.UNSCALED, method=PileupMethod.BAMNADO),
    log: OUTPUT_DIR + "/logs/heatmap/bamnado/unscaled/metaplot.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/bamnado/unscaled/metaplot.tsv",
    message: "Generating bamnado unscaled metaplot from matrix"


# ─── bamnado · merged · unscaled ─────────────────────────────────────────────

use rule heatmap_matrix as heatmap_matrix_bamnado_merged with:
    input:
        bigwigs=OUTPUT.select_bigwig_subtype(
            method=PileupMethod.BAMNADO,
            scale=DataScalingTechnique.UNSCALED,
            is_merged=True,
        ),
    output:
        matrix=temp(OUTPUT.select_heatmap_matrix(DataScalingTechnique.UNSCALED, method=PileupMethod.BAMNADO, is_merged=True)),
    log: OUTPUT_DIR + "/logs/heatmap/merged/bamnado/unscaled/matrix.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/merged/bamnado/unscaled/matrix.tsv",
    message: "Computing bamnado merged unscaled heatmap matrix from bigWig files"

use rule heatmap_plot as heatmap_plot_bamnado_merged with:
    input:
        matrix=OUTPUT.select_heatmap_matrix(DataScalingTechnique.UNSCALED, method=PileupMethod.BAMNADO, is_merged=True),
    output:
        heatmap=OUTPUT.select_heatmap_plot(DataScalingTechnique.UNSCALED, method=PileupMethod.BAMNADO, is_merged=True),
    log: OUTPUT_DIR + "/logs/heatmap/merged/bamnado/unscaled/heatmap.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/merged/bamnado/unscaled/heatmap.tsv",
    message: "Generating bamnado merged unscaled heatmap from matrix"

use rule heatmap_metaplot as heatmap_metaplot_bamnado_merged with:
    input:
        matrix=OUTPUT.select_heatmap_matrix(DataScalingTechnique.UNSCALED, method=PileupMethod.BAMNADO, is_merged=True),
    output:
        metaplot=OUTPUT.select_heatmap_metaplot(DataScalingTechnique.UNSCALED, method=PileupMethod.BAMNADO, is_merged=True),
    log: OUTPUT_DIR + "/logs/heatmap/merged/bamnado/unscaled/metaplot.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/merged/bamnado/unscaled/metaplot.tsv",
    message: "Generating bamnado merged unscaled metaplot from matrix"


# ─── bamnado · merged · csaw ─────────────────────────────────────────────────

use rule heatmap_matrix as heatmap_matrix_bamnado_merged_csaw with:
    input:
        bigwigs=OUTPUT.select_bigwig_subtype(
            method=PileupMethod.BAMNADO,
            scale=DataScalingTechnique.CSAW,
            is_merged=True,
        ),
    output:
        matrix=temp(OUTPUT.select_heatmap_matrix(DataScalingTechnique.CSAW, method=PileupMethod.BAMNADO, is_merged=True)),
    log: OUTPUT_DIR + "/logs/heatmap/merged/bamnado/csaw/matrix.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/merged/bamnado/csaw/matrix.tsv",
    message: "Computing bamnado merged CSAW-scaled heatmap matrix from bigWig files"

use rule heatmap_plot as heatmap_plot_bamnado_merged_csaw with:
    input:
        matrix=OUTPUT.select_heatmap_matrix(DataScalingTechnique.CSAW, method=PileupMethod.BAMNADO, is_merged=True),
    output:
        heatmap=OUTPUT.select_heatmap_plot(DataScalingTechnique.CSAW, method=PileupMethod.BAMNADO, is_merged=True),
    log: OUTPUT_DIR + "/logs/heatmap/merged/bamnado/csaw/heatmap.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/merged/bamnado/csaw/heatmap.tsv",
    message: "Generating bamnado merged CSAW-scaled heatmap from matrix"

use rule heatmap_metaplot as heatmap_metaplot_bamnado_merged_csaw with:
    input:
        matrix=OUTPUT.select_heatmap_matrix(DataScalingTechnique.CSAW, method=PileupMethod.BAMNADO, is_merged=True),
    output:
        metaplot=OUTPUT.select_heatmap_metaplot(DataScalingTechnique.CSAW, method=PileupMethod.BAMNADO, is_merged=True),
    log: OUTPUT_DIR + "/logs/heatmap/merged/bamnado/csaw/metaplot.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/merged/bamnado/csaw/metaplot.tsv",
    message: "Generating bamnado merged CSAW-scaled metaplot from matrix"


# ─── bamnado · merged · spikein ──────────────────────────────────────────────

use rule heatmap_matrix as heatmap_matrix_bamnado_merged_spikein with:
    input:
        bigwigs=OUTPUT.select_bigwig_subtype(
            method=PileupMethod.BAMNADO,
            scale=DataScalingTechnique.SPIKEIN,
            is_merged=True,
        ),
    output:
        matrix=temp(OUTPUT.select_heatmap_matrix(DataScalingTechnique.SPIKEIN, method=PileupMethod.BAMNADO, is_merged=True)),
    log: OUTPUT_DIR + "/logs/heatmap/merged/bamnado/spikein/matrix.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/merged/bamnado/spikein/matrix.tsv",
    message: "Computing bamnado merged spike-in normalized heatmap matrix from bigWig files"

use rule heatmap_plot as heatmap_plot_bamnado_merged_spikein with:
    input:
        matrix=OUTPUT.select_heatmap_matrix(DataScalingTechnique.SPIKEIN, method=PileupMethod.BAMNADO, is_merged=True),
    output:
        heatmap=OUTPUT.select_heatmap_plot(DataScalingTechnique.SPIKEIN, method=PileupMethod.BAMNADO, is_merged=True),
    log: OUTPUT_DIR + "/logs/heatmap/merged/bamnado/spikein/heatmap.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/merged/bamnado/spikein/heatmap.tsv",
    message: "Generating bamnado merged spike-in normalized heatmap from matrix"

use rule heatmap_metaplot as heatmap_metaplot_bamnado_merged_spikein with:
    input:
        matrix=OUTPUT.select_heatmap_matrix(DataScalingTechnique.SPIKEIN, method=PileupMethod.BAMNADO, is_merged=True),
    output:
        metaplot=OUTPUT.select_heatmap_metaplot(DataScalingTechnique.SPIKEIN, method=PileupMethod.BAMNADO, is_merged=True),
    log: OUTPUT_DIR + "/logs/heatmap/merged/bamnado/spikein/metaplot.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/merged/bamnado/spikein/metaplot.tsv",
    message: "Generating bamnado merged spike-in normalized metaplot from matrix"


# ─── homer · individual · unscaled ───────────────────────────────────────────

use rule heatmap_matrix as heatmap_matrix_homer with:
    input:
        bigwigs=OUTPUT.select_bigwig_subtype(
            method=PileupMethod.HOMER,
            scale=DataScalingTechnique.UNSCALED,
            ip_only=True,
        ),
    output:
        matrix=temp(OUTPUT.select_heatmap_matrix(DataScalingTechnique.UNSCALED, method=PileupMethod.HOMER)),
    log: OUTPUT_DIR + "/logs/heatmap/homer/unscaled/matrix.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/homer/unscaled/matrix.tsv",
    message: "Computing homer unscaled heatmap matrix from bigWig files"

use rule heatmap_plot as heatmap_plot_homer with:
    input:
        matrix=OUTPUT.select_heatmap_matrix(DataScalingTechnique.UNSCALED, method=PileupMethod.HOMER),
    output:
        heatmap=OUTPUT.select_heatmap_plot(DataScalingTechnique.UNSCALED, method=PileupMethod.HOMER),
    log: OUTPUT_DIR + "/logs/heatmap/homer/unscaled/heatmap.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/homer/unscaled/heatmap.tsv",
    message: "Generating homer unscaled heatmap from matrix"

use rule heatmap_metaplot as heatmap_metaplot_homer with:
    input:
        matrix=OUTPUT.select_heatmap_matrix(DataScalingTechnique.UNSCALED, method=PileupMethod.HOMER),
    output:
        metaplot=OUTPUT.select_heatmap_metaplot(DataScalingTechnique.UNSCALED, method=PileupMethod.HOMER),
    log: OUTPUT_DIR + "/logs/heatmap/homer/unscaled/metaplot.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/homer/unscaled/metaplot.tsv",
    message: "Generating homer unscaled metaplot from matrix"


# ─── homer · merged · unscaled ───────────────────────────────────────────────

use rule heatmap_matrix as heatmap_matrix_homer_merged with:
    input:
        bigwigs=OUTPUT.select_bigwig_subtype(
            method=PileupMethod.HOMER,
            scale=DataScalingTechnique.UNSCALED,
            is_merged=True,
        ),
    output:
        matrix=temp(OUTPUT.select_heatmap_matrix(DataScalingTechnique.UNSCALED, method=PileupMethod.HOMER, is_merged=True)),
    log: OUTPUT_DIR + "/logs/heatmap/merged/homer/unscaled/matrix.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/merged/homer/unscaled/matrix.tsv",
    message: "Computing homer merged unscaled heatmap matrix from bigWig files"

use rule heatmap_plot as heatmap_plot_homer_merged with:
    input:
        matrix=OUTPUT.select_heatmap_matrix(DataScalingTechnique.UNSCALED, method=PileupMethod.HOMER, is_merged=True),
    output:
        heatmap=OUTPUT.select_heatmap_plot(DataScalingTechnique.UNSCALED, method=PileupMethod.HOMER, is_merged=True),
    log: OUTPUT_DIR + "/logs/heatmap/merged/homer/unscaled/heatmap.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/merged/homer/unscaled/heatmap.tsv",
    message: "Generating homer merged unscaled heatmap from matrix"

use rule heatmap_metaplot as heatmap_metaplot_homer_merged with:
    input:
        matrix=OUTPUT.select_heatmap_matrix(DataScalingTechnique.UNSCALED, method=PileupMethod.HOMER, is_merged=True),
    output:
        metaplot=OUTPUT.select_heatmap_metaplot(DataScalingTechnique.UNSCALED, method=PileupMethod.HOMER, is_merged=True),
    log: OUTPUT_DIR + "/logs/heatmap/merged/homer/unscaled/metaplot.log",
    benchmark: OUTPUT_DIR + "/.benchmark/heatmap/merged/homer/unscaled/metaplot.tsv",
    message: "Generating homer merged unscaled metaplot from matrix"
