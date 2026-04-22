from seqnado.workflow.helpers.common import define_time_requested, define_memory_requested, get_read_count_flags

if CONFIG.qc.remove_blacklist:

    rule bam_remove_blacklisted_regions:
        input:
            bam=OUTPUT_DIR + "/aligned/sorted/{sample}.bam",
            bai=rules.bam_index.output.bai,
        output:
            bam=temp(OUTPUT_DIR + "/aligned/blacklist_regions_removed/{sample}.bam"),
            bai=temp(OUTPUT_DIR + "/aligned/blacklist_regions_removed/{sample}.bam.bai"),
        threads: 1
        params:
            blacklist=CONFIG.genome.blacklist,
            read_log=read_log_shared_path(OUTPUT_DIR, "{sample}"),
            count_flags=lambda wildcards: get_read_count_flags(wildcards, INPUT_FILES),
        resources:
            mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
            runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
        log: OUTPUT_DIR + "/logs/alignment_post_process/{sample}_blacklist.log",
        benchmark: OUTPUT_DIR + "/.benchmark/alignment_post_process/{sample}_blacklist.tsv",
        message: "Removing blacklisted regions from aligned BAM for sample {wildcards.sample} using bedtools",
        shell: f"""
        echo "Removing blacklisted regions with bedtools" > {{log}} 2>&1 &&
        before=$(samtools view -c {{params.count_flags}} {{input.bam}}) &&
        bedtools intersect -v -b {{params.blacklist}} -a {{input.bam}} > {{output.bam}} 2>> {{log}} &&
        samtools index -b {{output.bam}} -o {{output.bai}} >> {{log}} 2>&1 &&
        after=$(samtools view -c {{params.count_flags}} {{output.bam}}) &&
        {emit_read_logs("Blacklist", "{wildcards.sample}", "{params.read_log}")}
        """

else:

    rule bam_ignore_blacklisted_regions:
        input:
            bam=OUTPUT_DIR + "/aligned/sorted/{sample}.bam",
            bai=rules.bam_index.output.bai,
        output:
            bam=temp(OUTPUT_DIR + "/aligned/blacklist_regions_removed/{sample}.bam"),
            bai=temp(OUTPUT_DIR + "/aligned/blacklist_regions_removed/{sample}.bam.bai"),
        params:
            read_log=read_log_shared_path(OUTPUT_DIR, "{sample}"),
            count_flags=lambda wildcards: get_read_count_flags(wildcards, INPUT_FILES),
        threads: 1
        resources:
            mem=lambda wildcards, attempt: define_memory_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
            runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
        log: OUTPUT_DIR + "/logs/alignment_post_process/{sample}_blacklist.log",
        benchmark: OUTPUT_DIR + "/.benchmark/alignment_post_process/{sample}_blacklist.tsv",
        message: "Skipping blacklisted regions removal for sample {wildcards.sample}",
        shell: f"""
        echo "Skipping blacklisted regions removal" > {{log}} 2>&1 &&
        before=$(samtools view -c {{params.count_flags}} {{input.bam}}) &&
        mkdir -p $(dirname {{output.bam}}) &&
        cp {{input.bam}} {{output.bam}} >> {{log}} 2>&1 &&
        cp {{input.bai}} {{output.bai}} >> {{log}} 2>&1 &&
        after=$(samtools view -c {{params.count_flags}} {{output.bam}}) &&
        {emit_read_logs("Blacklist", "{wildcards.sample}", "{params.read_log}")}
        """
