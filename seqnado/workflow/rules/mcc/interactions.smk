import pandas as pd

from seqnado.workflow.helpers.mcc import viewpoint_to_grouped_viewpoint, extract_viewpoints
from seqnado.inputs.validation import ViewpointsFile

VIEWPOINTS_FILE = ViewpointsFile.validate(
    pd.read_csv(CONFIG.assay_config.mcc.viewpoints, 
    sep="\t", 
    names=["Chromosome", "Start", "End", "Name"]).assign(Score=0)
)
VIEWPOINT_OLIGOS = extract_viewpoints(CONFIG.assay_config.mcc.viewpoints)
VIEWPOINT_TO_GROUPED_VIEWPOINT = viewpoint_to_grouped_viewpoint(VIEWPOINT_OLIGOS)
GROUPED_VIEWPOINT_OLIGOS = list(set(VIEWPOINT_TO_GROUPED_VIEWPOINT.values()))

use rule sort_bam_by_qname as sort_mcc_annotated_bam with:
    input:
        bam=OUTPUT_DIR + "/mcc/{group}/{group}.bam",
    output:
        bam=temp(OUTPUT_DIR + "/mcc/{group}/{group}_qname.bam"),
        read_log=temp(OUTPUT_DIR + "/qc/mcc/{group}_qname_sort.tsv"),
    threads: CONFIG.third_party_tools.samtools.sort.threads
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/mcc/sort_genomic_aligned_reads/{group}.log",
    message: "Sorting aligned BAM by QNAME for group {wildcards.group} using samtools",
    benchmark: OUTPUT_DIR + "/.benchmark/mcc/sort_genomic_aligned_reads/{group}.tsv",


rule identify_ligation_junctions_grouped:
    input:
        bam=OUTPUT_DIR + "/mcc/{group}/{group}_qname.bam"
    output:
        pairs=expand(OUTPUT_DIR + "/mcc/{{group}}/ligation_junctions/raw/{viewpoint}.pairs", viewpoint=GROUPED_VIEWPOINT_OLIGOS),
    params:
        outdir=OUTPUT_DIR + "/mcc/{group}/ligation_junctions/raw/",
    threads: 1
    resources:
        mem="1GB",
    container: 'docker://ghcr.io/alsmith151/mccnado:latest',
    wildcard_constraints:
        group="|".join(SAMPLE_GROUPINGS.get_grouping('consensus').group_names),
    log: OUTPUT_DIR + "/logs/ligation_junctions/{group}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/ligation_junctions/{group}.tsv",
    message: "Identifying ligation junctions for group {wildcards.group}",
    shell: """
    mccnado identify-ligation-junctions \
    {input.bam} \
    {params.outdir}
    """


rule sort_ligation_junctions:
    input:
        pairs=OUTPUT_DIR + "/mcc/{group}/ligation_junctions/raw/{viewpoint}.pairs",
    output:
        pairs=temp(OUTPUT_DIR + "/mcc/{group}/ligation_junctions/{viewpoint}.pairs"),
    threads: 1
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
    log: OUTPUT_DIR + "/logs/sort_ligation_junctions/{group}_{viewpoint}.log",
    wildcard_constraints:
        group="|".join(SAMPLE_GROUPINGS.get_grouping('consensus').group_names),
    benchmark: OUTPUT_DIR + "/.benchmark/sort_ligation_junctions/{group}_{viewpoint}.tsv",
    message: "Sorting ligation junctions for viewpoint {wildcards.viewpoint} in group {wildcards.group}",
    shell: """
    sort -k2,2 -k4,4 -k3,3n -k5,5n {input.pairs} > {output.pairs}
    """


rule bgzip_pairs:
    input:
        pairs=OUTPUT_DIR + "/mcc/{group}/ligation_junctions/{viewpoint}.pairs",
    output:
        pairs=OUTPUT_DIR + "/mcc/{group}/ligation_junctions/{viewpoint}.pairs.gz",
    log: OUTPUT_DIR + "/logs/bgzip_pairs/{group}_{viewpoint}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/bgzip_pairs/{group}_{viewpoint}.tsv",
    message: "Bgzipping pairs file for viewpoint {wildcards.viewpoint} in group {wildcards.group}",
    wildcard_constraints:
        group="|".join(SAMPLE_GROUPINGS.get_grouping('consensus').group_names),
    shell: """
    bgzip -c {input.pairs} > {output.pairs}
    """


rule identify_ligation_junctions_replicates:
    input:
        bam=OUTPUT_DIR + "/mcc/replicates/{sample}/{sample}_deduplicated.bam",
    output:
        pairs=expand(OUTPUT_DIR + "/mcc/replicates/{{sample}}/ligation_junctions/{viewpoint}.pairs", viewpoint=GROUPED_VIEWPOINT_OLIGOS),
    params:
        outdir=OUTPUT_DIR + "/mcc/replicates/{sample}/ligation_junctions/",
    threads: 1
    resources:
        mem="1GB",
    container: 'docker://ghcr.io/alsmith151/mccnado:latest',
    log: OUTPUT_DIR + "/logs/ligation_junctions/{sample}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/ligation_junctions/{sample}.tsv",
    message: "Identifying ligation junctions for sample {wildcards.sample}",
    shell: """
    mccnado identify-ligation-junctions \
    {input.bam} \
    {params.outdir}
    """

use rule sort_ligation_junctions as sort_ligation_junctions_replicates with:
    input:
        pairs=OUTPUT_DIR + "/mcc/replicates/{sample}/ligation_junctions/{viewpoint}.pairs",
    output:
        pairs=temp(OUTPUT_DIR + "/mcc/replicates/{sample}/ligation_junctions/sorted/{viewpoint}.pairs"),
    threads: 1
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
    log: OUTPUT_DIR + "/logs/sort_ligation_junctions/{sample}_{viewpoint}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/sort_ligation_junctions/{sample}_{viewpoint}.tsv",
    message: "Sorting ligation junctions for viewpoint {wildcards.viewpoint} in sample {wildcards.sample}",

use rule bgzip_pairs as bgzip_pairs_replicates with:
    input:
        pairs=OUTPUT_DIR + "/mcc/replicates/{sample}/ligation_junctions/sorted/{viewpoint}.pairs",
    output:
        pairs=OUTPUT_DIR + "/mcc/replicates/{sample}/ligation_junctions/{viewpoint}.pairs.gz",
    log: OUTPUT_DIR + "/logs/bgzip_pairs/{sample}_{viewpoint}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/bgzip_pairs/{sample}_{viewpoint}.tsv",
    message: "Bgzipping pairs file for viewpoint {wildcards.viewpoint} in sample {wildcards.sample}",


rule make_cooler:
    input:
        pairs=OUTPUT_DIR + "/mcc/{group}/ligation_junctions/{viewpoint}.pairs.gz",
        bins=OUTPUT_DIR + "/resources/genomic_bins.bed",
    output:
        cooler=temp(OUTPUT_DIR + "/mcc/{group}/ligation_junctions/{viewpoint}.cool"),
    params:
        resolution=CONFIG.assay_config.mcc.resolutions[0],
        genome=CONFIG.genome.name,
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    wildcard_constraints:
        group="|".join(SAMPLE_GROUPINGS.get_grouping('consensus').group_names),
    resources:
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=8, attempts=attempt, scale=SCALE_RESOURCES),
    log: OUTPUT_DIR + "/logs/make_cooler/{group}_{viewpoint}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/make_cooler/{group}_{viewpoint}.tsv",
    message: "Creating cooler file for viewpoint {wildcards.viewpoint} in group {wildcards.group}",
    shell: """
    cooler cload pairs \
    {input.bins} \
    {input.pairs} \
    {output.cooler} \
    --assembly {params.genome} \
    -c1 2 -p1 3 -c2 4 -p2 5 > {log} 2>&1
    """


rule zoomify_cooler:
    input:
        cooler=OUTPUT_DIR + "/mcc/{group}/ligation_junctions/{viewpoint}.cool",
    output:
        cooler=temp(OUTPUT_DIR + "/mcc/{group}/ligation_junctions/{viewpoint}.mcool"),
    params:
        resolutions=[f'-r {res}' for res in CONFIG.assay_config.mcc.resolutions],
    resources:
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=8, attempts=attempt, scale=SCALE_RESOURCES),
    wildcard_constraints:
        group="|".join(SAMPLE_GROUPINGS.get_grouping('consensus').group_names),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/zoomify_cooler/{group}_{viewpoint}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/zoomify_cooler/{group}_{viewpoint}.tsv",
    message: "Zoomifying cooler file for viewpoint {wildcards.viewpoint} in group {wildcards.group}",
    shell: """
    cooler zoomify {input.cooler} {params.resolutions} -o {output.cooler} > {log} 2>&1
    """


rule aggregate_coolers:
    input:
        mcools=expand(OUTPUT_DIR + "/mcc/{{group}}/ligation_junctions/{viewpoint}.mcool", 
                    viewpoint=GROUPED_VIEWPOINT_OLIGOS),
    output:
        mcool=OUTPUT_DIR + "/mcc/contacts/{group}/{group}.mcool",
    resources:
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=8, attempts=attempt, scale=SCALE_RESOURCES),
    container: "docker://ghcr.io/alsmith151/mccnado:latest"
    log: OUTPUT_DIR + "/logs/{group}_aggregate_coolers.log",
    benchmark: OUTPUT_DIR + "/.benchmark/{group}_aggregate_coolers.tsv",
    message: "Aggregating cooler files for group {wildcards.group}",
    shell: """
        mccnado combine-ligation-junction-coolers \
        {input.mcools} \
        {output.mcool}
        """


rule create_sentinel_contact_files:
    input:
        mcool_grouped=expand(OUTPUT_DIR + "/mcc/contacts/{group}/{group}.mcool", 
                    group=SAMPLE_GROUPINGS.get_grouping('consensus').group_names),
        mcool_replicate=expand(OUTPUT_DIR + "/mcc/replicates/{sample}/ligation_junctions/{viewpoint}.pairs.gz", 
                    viewpoint=GROUPED_VIEWPOINT_OLIGOS,
                    sample=SAMPLE_NAMES) if CONFIG.assay_config.mcc.create_replicate_contact_files else [],
    output:
        sentinel_contacts=OUTPUT_DIR + "/mcc/.mcc_contacts_identified.txt",
    resources:
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=8, attempts=attempt, scale=SCALE_RESOURCES),
    shell: """
        touch {output.sentinel_contacts}
        """