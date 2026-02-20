from seqnado import Assay

rule save_design:
    output:
        OUTPUT_DIR + "/metadata.csv",
    container: None
    log: OUTPUT_DIR + "/logs/save_design.log",
    benchmark: OUTPUT_DIR + "/.benchmark/save_design.tsv",
    message: "Saving design dataframe to metadata.csv in output directory",
    run:
        INPUT_FILES.to_dataframe().to_csv(OUTPUT_DIR + "/metadata.csv", index=False)


rule make_genomic_bins:
    input:
        chrom_sizes=CONFIG.genome.chromosome_sizes,
    params:
        bin_size=CONFIG.genome.bin_size if not ASSAY == Assay.MCC else CONFIG.assay_config.mcc.resolutions[0],
    output:
        bed=OUTPUT_DIR + "/resources/genomic_bins.bed",
    container:
        "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    resources:
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    log: OUTPUT_DIR + "/logs/genomic_bins.log",
    benchmark: OUTPUT_DIR + "/.benchmark/genomic_bins.tsv",
    message: "Generating genomic bins of size {params.bin_size} from chromosome sizes",
    shell: """
    cooler makebins {input.chrom_sizes} {params.bin_size} -o {output.bed} > {log} 2>&1
    """


rule bed_to_saf:
    input:
        bed=OUTPUT_DIR + "/resources/genomic_bins.bed",
    output:
        saf=OUTPUT_DIR + "/resources/genomic_bins.saf",
    container:
        "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    resources:
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    log:
        OUTPUT_DIR + "/logs/genomic_bins_saf.log",
    benchmark: OUTPUT_DIR + "/.benchmark/genomic_bins_saf.tsv",
    message: "Converting genomic bins BED to SAF format",
    shell: """
    awk 'BEGIN{{OFS="\t"}} {{print $1":"$2"-"$3, $1, $2, $3, "."}}' {input.bed} > {output.saf} 2> {log}
    """


rule validate_peaks:
    input:
        peaks=OUTPUT.peak_files,
    output:
        sentinel=OUTPUT_DIR + "/peaks/.validated",
    container:
        "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log:
        OUTPUT_DIR + "/logs/validate_peaks.log",
    benchmark: OUTPUT_DIR + "/.benchmark/validate_peaks.tsv",
    message: "Validating peak files to ensure they are not empty",
    run:
        from loguru import logger

        with logger.catch():
            for peak_file in input.peaks:
                with open(peak_file, "r+") as p:
                    peak_entries = p.readlines()
                    n_peak_lines = sum(1 for line in peak_entries if not line.startswith("#"))
                    if n_peak_lines == 0: # empty peak file, write a dummy peak
                        p.write("chr21\t1\t2\n")

        with open(output.sentinel, "w") as s:
            s.write("validated")


rule bed_to_bigbed:
    input:
        bed=OUTPUT_DIR + "/peaks/{directory}/{sample}.bed",
        sentinel=OUTPUT_DIR + "/peaks/.validated",
    output:
        bigbed=OUTPUT_DIR + "/peaks/{directory}/{sample}.bigBed",
    params:
        chrom_sizes=CONFIG.genome.chromosome_sizes,
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/bed_to_bigbed/{directory}/{sample}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/bed_to_bigbed/{directory}/{sample}.tsv",
    message: "Converting BED to BigBed for sample {wildcards.sample} in directory {wildcards.directory}",
    shell:
        """
        sort -k1,1 -k2,2n {input.bed} | grep '#' -v | cut -f 1-4 > {input.bed}.tmp &&
        bedToBigBed {input.bed}.tmp {params.chrom_sizes} {output.bigbed} 2> {log} &&
        rm {input.bed}.tmp
        """



localrules:
    save_design,
    make_genomic_bins,
    bed_to_saf,
    validate_peaks,
