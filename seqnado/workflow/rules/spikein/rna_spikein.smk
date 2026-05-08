from seqnado.workflow.helpers.common import define_memory_requested, define_time_requested

use rule align_paired as align_paired_spikein_rna with:
    output:
        bam_dir=temp(directory(OUTPUT_DIR + "/aligned/spikein/{sample}.star/")),
        log_out=temp(OUTPUT_DIR + "/aligned/spikein/{sample}.star/Log.final.out"),
    params:
        index=str(CONFIG.genome.index.prefix),
        options=str(CONFIG.third_party_tools.star.align.command_line_arguments),
        prefix=OUTPUT_DIR + "/aligned/spikein/{sample}.star/",
        read_log=read_log_shared_path(OUTPUT_DIR, "{sample}"),
        rule_label="align_paired_spikein_rna",
    log: OUTPUT_DIR + "/logs/align_spikein/{sample}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/align_spikein/{sample}.tsv",

use rule align_single as align_single_spikein_rna with:
    output:
        bam_dir=temp(directory(OUTPUT_DIR + "/aligned/spikein/{sample}.star/")),
        log_out=temp(OUTPUT_DIR + "/aligned/spikein/{sample}.star/Log.final.out"),
    params:
        index=str(CONFIG.genome.index.prefix),
        options=str(CONFIG.third_party_tools.star.align.command_line_arguments),
        prefix=OUTPUT_DIR + "/aligned/spikein/{sample}.star/",
        read_log=read_log_shared_path(OUTPUT_DIR, "{sample}"),
        rule_label="align_single_spikein_rna",
    log: OUTPUT_DIR + "/logs/align_spikein/{sample}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/align_spikein/{sample}.tsv",

use rule rename_aligned as rename_aligned_spikein with:
    input:
        bam_dir=OUTPUT_DIR + "/aligned/spikein/{sample}.star/",
        log_out=OUTPUT_DIR + "/aligned/spikein/{sample}.star/Log.final.out",
    output:
        bam=temp(OUTPUT_DIR + "/aligned/spikein/raw/{sample}.bam"),
    params:
        read_log=read_log_shared_path(OUTPUT_DIR, "{sample}"),
