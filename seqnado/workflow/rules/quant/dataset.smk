from seqnado.workflow.helpers.common import define_time_requested, define_memory_requested


rule make_dataset_regions:
    """Create a dataset from bigWig files using either a BED file."""
    input:
        bigwigs=OUTPUT.select_bigwig_subtype(
            method=PileupMethod.DEEPTOOLS,
            scale=ScaleMethod.UNSCALED
        ),
    output:
        dataset=OUTPUT_DIR + "/dataset/dataset_regions.h5ad",
    params:
        bigwig_dir=OUTPUT_DIR + "/bigwigs/deeptools/",
        chromosome_sizes=CONFIG.genome.chromosome_sizes,
        blacklist=CONFIG.genome.blacklist,
        regions=CONFIG.assay_config.dataset_for_ml.regions_bed,
    threads: 1
    resources:
            mem=lambda wildcards, attempt: define_memory_requested(initial_value=32, attempts=attempt, scale=SCALE_RESOURCES),
            runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_ml_cpu:latest",
    log: OUTPUT_DIR + "/logs/make_dataset_regions.log",
    benchmark: OUTPUT_DIR + "/.benchmark/make_dataset_regions.tsv",
    message: "Making dataset from regions for machine learning"
    shell: """
    quantnado-make-dataset \
    --bigwig-dir {params.bigwig_dir} \
    --output-file {output.dataset} \
    --chromsizes {params.chromosome_sizes} \
    --regions {params.regions} \
    --blacklist {params.blacklist} \
    --log-file {log}
    """


rule make_dataset_binsize:
    """Create a dataset from bigWig files using bin size."""
    input:
        bigwigs=OUTPUT.select_bigwig_subtype(
            method=PileupMethod.DEEPTOOLS,
            scale=ScaleMethod.UNSCALED
        ),
    output:
        dataset=OUTPUT_DIR + "/dataset/dataset_bins.h5ad",
    params:
        bigwig_dir=OUTPUT_DIR + "/bigwigs/deeptools/",
        chromosome_sizes=CONFIG.genome.chromosome_sizes,
        blacklist=CONFIG.genome.blacklist,
        binsize=CONFIG.assay_config.dataset_for_ml.binsize,
    threads: 1
    resources:
            mem=lambda wildcards, attempt: define_memory_requested(initial_value=32, attempts=attempt, scale=SCALE_RESOURCES),
            runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_ml_cpu:latest",
    log: OUTPUT_DIR + "/logs/make_dataset_binsize.log",
    benchmark: OUTPUT_DIR + "/.benchmark/make_dataset_binsize.tsv",
    message: "Making dataset from binsize for machine learning"
    shell: """
    quantnado-make-dataset \
    --bigwig-dir {params.bigwig_dir} \
    --output-file {output.dataset} \
    --chromsizes {params.chromosome_sizes} \
    --binsize {params.binsize} \
    --blacklist {params.blacklist} \
    --log-file {log}
    """
