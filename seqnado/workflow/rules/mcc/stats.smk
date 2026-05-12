rule extract_ligation_stats:
    input:
        bam=OUTPUT_DIR + "/mcc/replicates/{sample}/{sample}.bam",
    output:
        stats=OUTPUT_DIR + "/resources/replicates/{sample}_ligation_stats.json"
    container: 'docker://ghcr.io/alsmith151/mccnado:latest'
    log: OUTPUT_DIR + "/logs/extract_ligation_stats/{sample}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/extract_ligation_stats/{sample}.tsv",
    message: "Extracting ligation stats for MCC BAM of sample {wildcards.sample}",
    shell: """
    mccnado extract-ligation-stats {input.bam} {output.stats} 
    """

use rule extract_ligation_stats as extract_ligation_stats_merged with:
    input:
        bam=OUTPUT_DIR + "/mcc/{group}/{group}.bam",
    output:
        stats=OUTPUT_DIR + "/resources/{group}_ligation_stats.json",
    wildcard_constraints:
        group="|".join(SAMPLE_GROUPINGS.get_grouping('consensus').group_names),
    log: OUTPUT_DIR + "/logs/extract_ligation_stats_merged/{group}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/extract_ligation_stats_merged/{group}.tsv",
    message: "Extracting ligation stats for merged MCC BAM of group {wildcards.group}",
    