from seqnado import ScalingMethod
from seqnado.workflow.helpers.bam import get_bam_files_for_scaling_group
from seqnado.workflow.helpers.common import define_time_requested, define_memory_requested


rule bamnado_scaling_factors:
    input:
        bam=lambda wc: get_bam_files_for_scaling_group(wc, SAMPLE_GROUPINGS, OUTPUT_DIR),
        bai=lambda wc: [b + ".bai" for b in get_bam_files_for_scaling_group(wc, SAMPLE_GROUPINGS, OUTPUT_DIR)],
    output:
        scaling_factors=OUTPUT_DIR + "/resources/{scaling_method}/{group}_scaling_factors.tsv",
    params:
        bams=lambda wc, input: " ".join(f"-b {b}" for b in input.bam),
        method=lambda wc: wc.scaling_method.replace("_", "-"),
        options=str(CONFIG.third_party_tools.bamnado.bam_normalize.command_line_arguments),
    threads: CONFIG.third_party_tools.bamnado.bam_normalize.threads
    wildcard_constraints:
        scaling_method="|".join(m.value for m in ScalingMethod),
        group="|".join(
            {
                *SAMPLE_GROUPINGS.get_grouping("scaling").group_names,
                *SAMPLE_GROUPINGS.get_grouping("consensus").group_names,
            }
        ),
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
    container: "docker://ghcr.io/alsmith151/bamnado:latest"
    log:
        OUTPUT_DIR + "/logs/normalization/{scaling_method}/{group}_scaling_factors.log"
    benchmark:
        OUTPUT_DIR + "/.benchmark/normalization/{scaling_method}/{group}_scaling_factors.tsv"
    message:
        "Calculating {wildcards.scaling_method} scaling factors for group {wildcards.group}"
    shell: """
    bamnado bam-normalize {params.bams} --method {params.method} --format tsv {params.options} -o {output.scaling_factors} > {log} 2>&1
    """
