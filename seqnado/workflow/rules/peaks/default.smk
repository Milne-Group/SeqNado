from seqnado import FileType
from seqnado.workflow.helpers.peaks import (
    get_control_file,
    correct_macs_options,
    get_lanceotron_call_peaks_threshold,
)


rule call_peaks_macs2_with_input:
    input:
        treatment=OUTPUT_DIR + "/aligned/{sample_id}.bam",
        control=lambda wc: get_control_file(wc, file_type=FileType.BAM, INPUT_FILES=INPUT_FILES, OUTPUT_DIR=OUTPUT_DIR),
    output:
        peaks=OUTPUT_DIR + "/peaks/macs2/{sample_id}.bed",
    wildcard_constraints:
        sample_id=r"(?!merged/).*",
    params:
        options=lambda wc: str(
            correct_macs_options(
                wc, CONFIG.third_party_tools.macs.call_peaks.command_line_arguments,
                INPUT_FILES=INPUT_FILES, SAMPLE_GROUPINGS=SAMPLE_GROUPINGS
            )
        ),
        raw=lambda wc, output: output.peaks.replace(".bed", "_peaks.xls"),
        basename=lambda wc, output: output.peaks.replace(".bed", ""),
        narrow_peak=lambda wc, output: output.peaks.replace(".bed", "_peaks.narrowPeak"),
    threads: 1
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(
            initial_value=4, attempts=attempt, scale=SCALE_RESOURCES
        ),
        runtime=lambda wildcards, attempt: define_time_requested(
            initial_value=6, attempts=attempt, scale=SCALE_RESOURCES
        ),
    container:
        "docker://quay.io/biocontainers/macs2:2.1.1.20160309--py27r3.3.1_1"
    log:
        OUTPUT_DIR + "/logs/macs2/{sample_id}.log",
    benchmark:
        OUTPUT_DIR + "/.benchmark/peaks/macs2/{sample_id}.tsv"
    message:
        "Calling peaks with MACS2 for sample {wildcards.sample_id}"
    shell:
        """
        if ! macs2 callpeak -t {input.treatment} -c {input.control} -n {params.basename} {params.options} > {log} 2>&1; then
            touch {output.peaks}
        else
            awk 'BEGIN{{OFS="\\t"}} !/^#/ && !/^chr[[:space:]]+start[[:space:]]+end/ && !/^$/ {{print $1, $2, $3}}' {params.narrow_peak} > {output.peaks} 2>> {log}
        fi
        """


rule call_peaks_macs2_no_input:
    input:
        treatment=OUTPUT_DIR + "/aligned/{sample_id}.bam",
    output:
        peaks=OUTPUT_DIR + "/peaks/macs2/{sample_id}.bed",
    wildcard_constraints:
        sample_id=r"(?!merged/).*",
    params:
        options=lambda wc: str(
            correct_macs_options(
                wc, CONFIG.third_party_tools.macs.call_peaks.command_line_arguments,
                INPUT_FILES=INPUT_FILES, SAMPLE_GROUPINGS=SAMPLE_GROUPINGS
            )
        ),
        raw=lambda wc, output: output.peaks.replace(".bed", "_peaks.xls"),
        basename=lambda wc, output: output.peaks.replace(".bed", ""),
        narrow_peak=lambda wc, output: output.peaks.replace(".bed", "_peaks.narrowPeak"),
    threads: 1
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(
            initial_value=2, attempts=attempt, scale=SCALE_RESOURCES
        ),
        runtime=lambda wildcards, attempt: define_time_requested(
            initial_value=2, attempts=attempt, scale=SCALE_RESOURCES
        ),
    container:
        "docker://quay.io/biocontainers/macs2:2.1.1.20160309--py27r3.3.1_1"
    log:
        OUTPUT_DIR + "/logs/macs2/{sample_id}.log",
    benchmark:
        OUTPUT_DIR + "/.benchmark/peaks/macs2/{sample_id}.tsv"
    message:
        "Calling peaks with MACS2 for sample {wildcards.sample_id}"
    shell:
        """
        if ! macs2 callpeak -t {input.treatment} -n {params.basename} {params.options} > {log} 2>&1; then
            touch {output.peaks}
        else
            awk 'BEGIN{{OFS="\\t"}} !/^#/ && !/^chr[[:space:]]+start[[:space:]]+end/ && !/^$/ {{print $1, $2, $3}}' {params.narrow_peak} > {output.peaks} 2>> {log}
        fi
        """


ruleorder: call_peaks_macs2_with_input > call_peaks_macs2_no_input


# MACS3 rules (same as MACS2 but with macs3 container)


rule call_peaks_macs3_with_input:
    input:
        treatment=OUTPUT_DIR + "/aligned/{sample_id}.bam",
        control=lambda wc: get_control_file(wc, file_type=FileType.BAM, INPUT_FILES=INPUT_FILES, OUTPUT_DIR=OUTPUT_DIR),
    output:
        peaks=OUTPUT_DIR + "/peaks/macs3/{sample_id}.bed",
    wildcard_constraints:
        sample_id=r"(?!merged/).*",
    params:
        options=lambda wc: str(
            correct_macs_options(
                wc, CONFIG.third_party_tools.macs.call_peaks.command_line_arguments,
                INPUT_FILES=INPUT_FILES, SAMPLE_GROUPINGS=SAMPLE_GROUPINGS
            )
        ),
        raw=lambda wc, output: output.peaks.replace(".bed", "_peaks.xls"),
        basename=lambda wc, output: output.peaks.replace(".bed", ""),
        narrow_peak=lambda wc, output: output.peaks.replace(".bed", "_peaks.narrowPeak"),
    threads: 1
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(
            initial_value=4, attempts=attempt, scale=SCALE_RESOURCES
        ),
        runtime=lambda wildcards, attempt: define_time_requested(
            initial_value=6, attempts=attempt, scale=SCALE_RESOURCES
        ),
    container:
        "docker://quay.io/biocontainers/macs3:3.0.3--py39h0699b22_0"
    log:
        OUTPUT_DIR + "/logs/macs3/{sample_id}.log",
    benchmark:
        OUTPUT_DIR + "/.benchmark/peaks/macs3/{sample_id}.tsv"
    message:
        "Calling peaks with MACS3 for sample {wildcards.sample_id}"
    shell:
        """
        if ! macs3 callpeak -t {input.treatment} -c {input.control} -n {params.basename} {params.options} > {log} 2>&1; then
            touch {output.peaks}
        else
            awk 'BEGIN{{OFS="\\t"}} !/^#/ && !/^chr[[:space:]]+start[[:space:]]+end/ && !/^$/ {{print $1, $2, $3}}' {params.narrow_peak} > {output.peaks} 2>> {log}
        fi
        """


rule call_peaks_macs3_no_input:
    input:
        treatment=OUTPUT_DIR + "/aligned/{sample_id}.bam",
    output:
        peaks=OUTPUT_DIR + "/peaks/macs3/{sample_id}.bed",
    wildcard_constraints:
        sample_id=r"(?!merged/).*",
    params:
        options=lambda wc: str(
            correct_macs_options(
                wc, CONFIG.third_party_tools.macs.call_peaks.command_line_arguments,
                INPUT_FILES=INPUT_FILES, SAMPLE_GROUPINGS=SAMPLE_GROUPINGS
            )
        ),
        raw=lambda wc, output: output.peaks.replace(".bed", "_peaks.xls"),
        basename=lambda wc, output: output.peaks.replace(".bed", ""),
        narrow_peak=lambda wc, output: output.peaks.replace(".bed", "_peaks.narrowPeak"),
    threads: 1
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(
            initial_value=2, attempts=attempt, scale=SCALE_RESOURCES
        ),
        runtime=lambda wildcards, attempt: define_time_requested(
            initial_value=2, attempts=attempt, scale=SCALE_RESOURCES
        ),
    container:
        "docker://quay.io/biocontainers/macs3:3.0.3--py39h0699b22_0"
    log:
        OUTPUT_DIR + "/logs/macs3/{sample_id}.log",
    benchmark:
        OUTPUT_DIR + "/.benchmark/peaks/macs3/{sample_id}.tsv"
    message:
        "Calling peaks with MACS3 for sample {wildcards.sample_id}"
    shell:
        """
        if ! macs3 callpeak -t {input.treatment} -n {params.basename} {params.options} > {log} 2>&1; then
            touch {output.peaks}
        else
            awk 'BEGIN{{OFS="\\t"}} !/^#/ && !/^chr[[:space:]]+start[[:space:]]+end/ && !/^$/ {{print $1, $2, $3}}' {params.narrow_peak} > {output.peaks} 2>> {log}
        fi
        """


ruleorder: call_peaks_macs3_with_input > call_peaks_macs3_no_input


rule call_peaks_homer_with_input:
    input:
        treatment=OUTPUT_DIR + "/tag_dirs/{sample_id}",
        control=lambda wc: get_control_file(wc, file_type=FileType.TAG_DIRECTORY, INPUT_FILES=INPUT_FILES, OUTPUT_DIR=OUTPUT_DIR),
    output:
        peaks=OUTPUT_DIR + "/peaks/homer/{sample_id}.bed",
    wildcard_constraints:
        sample_id=r"(?!merged/).*",
    params:
        options=str(CONFIG.third_party_tools.homer.find_peaks.command_line_arguments),
    threads: 1
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(
            initial_value=4, attempts=attempt, scale=SCALE_RESOURCES
        ),
        runtime=lambda wildcards, attempt: define_time_requested(
            initial_value=2, attempts=attempt, scale=SCALE_RESOURCES
        ),
    container:
        "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log:
        OUTPUT_DIR + "/logs/homer/{sample_id}.log",
    benchmark:
        OUTPUT_DIR + "/.benchmark/homer/{sample_id}.tsv"
    message:
        "Calling peaks with HOMER for sample {wildcards.sample_id}"
    shell:
        """
        findPeaks {input.treatment} {params.options} -o {output.peaks}.tmp  -i {input.control} > {log} 2>&1 &&
        pos2bed.pl {output.peaks}.tmp -o {output.peaks} >> {log} 2>&1 &&
        rm {output.peaks}.tmp
        """


rule call_peaks_homer_no_input:
    input:
        treatment=OUTPUT_DIR + "/tag_dirs/{sample_id}",
    output:
        peaks=OUTPUT_DIR + "/peaks/homer/{sample_id}.bed",
    wildcard_constraints:
        sample_id=r"(?!merged/).*",
    params:
        options=str(CONFIG.third_party_tools.homer.find_peaks.command_line_arguments),
    threads: 1
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(
            initial_value=4, attempts=attempt, scale=SCALE_RESOURCES
        ),
        runtime=lambda wildcards, attempt: define_time_requested(
            initial_value=2, attempts=attempt, scale=SCALE_RESOURCES
        ),
    container:
        "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log:
        OUTPUT_DIR + "/logs/homer/{sample_id}.log",
    benchmark:
        OUTPUT_DIR + "/.benchmark/homer/{sample_id}.tsv"
    message:
        "Calling peaks with HOMER for sample {wildcards.sample_id}"
    shell:
        """
        findPeaks {input.treatment} {params.options} -o {output.peaks}.tmp > {log} 2>&1 &&
        pos2bed.pl {output.peaks}.tmp -o {output.peaks} >> {log} 2>&1 &&
        rm {output.peaks}.tmp
        """


ruleorder: call_peaks_homer_with_input > call_peaks_homer_no_input


rule call_peaks_lanceotron_with_input:
    input:
        treatment=OUTPUT_DIR + "/bigwigs/deeptools/unscaled/{sample_id}.bigWig",
        control=lambda wc: get_control_file(wc, file_type=FileType.BIGWIG, INPUT_FILES=INPUT_FILES, OUTPUT_DIR=OUTPUT_DIR),
    output:
        peaks=OUTPUT_DIR + "/peaks/lanceotron/{sample_id}.bed",
        ltron_peaks=temp(OUTPUT_DIR + "/peaks/lanceotron/{sample_id}_L-tron.bed"),
    wildcard_constraints:
        sample_id=r"(?!merged/).*",
    params:
        threshold=lambda wildcards: get_lanceotron_call_peaks_threshold(wildcards, CONFIG=CONFIG),
        options=str(
            CONFIG.third_party_tools.lanceotron.call_peaks.command_line_arguments
        ),
        outdir=lambda wc, output: os.path.dirname(output.peaks),
        basename=lambda wc, output: output.peaks.replace(".bed", ""),
    container:
        "oras://ghcr.io/alsmith151/seqnado_ml_cpu:latest"
    threads: 
        CONFIG.third_party_tools.lanceotron.call_peaks.threads
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(
            initial_value=16, attempts=attempt, scale=SCALE_RESOURCES
        ),
        runtime=lambda wildcards, attempt: define_time_requested(
            initial_value=6, attempts=attempt, scale=SCALE_RESOURCES
        ),
    log:
        OUTPUT_DIR + "/logs/lanceotron/{sample_id}.log",
    benchmark:
        OUTPUT_DIR + "/.benchmark/lanceotron/{sample_id}.tsv"
    message:
        "Calling peaks with LanceOtron for sample {wildcards.sample_id}"
    shell:
        """
        lanceotron callPeaksInput {input.treatment} -i {input.control} -f {params.outdir} --skipheader > {log} 2>&1 &&
        cat {output.ltron_peaks} | awk 'BEGIN{{OFS="\\t"}} $4 >= {params.threshold} {{print $1, $2, $3}}' > {output.peaks} 
        """


rule call_peaks_lanceotron_no_input:
    input:
        treatment=OUTPUT_DIR + "/bigwigs/deeptools/unscaled/{sample_id}.bigWig",
    output:
        peaks=OUTPUT_DIR + "/peaks/lanceotron/{sample_id}.bed",
        ltron_peaks=temp(OUTPUT_DIR + "/peaks/lanceotron/{sample_id}_L-tron.bed"),
    wildcard_constraints:
        sample_id=r"(?!merged/).*",
    params:
        options=str(
            CONFIG.third_party_tools.lanceotron.call_peaks.command_line_arguments
        ),
        outdir=lambda wc, output: os.path.dirname(output.peaks),
        basename=lambda wc, output: output.peaks.replace(".bed", ""),
    threads: 
        CONFIG.third_party_tools.lanceotron.call_peaks.threads
    container:
        "oras://ghcr.io/alsmith151/seqnado_ml_cpu:latest"
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(
            initial_value=16, attempts=attempt, scale=SCALE_RESOURCES
        ),
        runtime=lambda wildcards, attempt: define_time_requested(
            initial_value=6, attempts=attempt, scale=SCALE_RESOURCES
        ),
    log:
        OUTPUT_DIR + "/logs/lanceotron/{sample_id}.log",
    benchmark:
        OUTPUT_DIR + "/.benchmark/lanceotron/{sample_id}.tsv"
    message:
        "Calling peaks with LanceOtron for sample {wildcards.sample_id}"
    shell:
        """
        lanceotron callPeaks {input.treatment} -f {params.outdir} --skipheader  {params.options} > {log} 2>&1 &&
        cat {output.ltron_peaks} | cut -f 1-3 > {output.peaks}
        """


ruleorder: call_peaks_lanceotron_with_input > call_peaks_lanceotron_no_input


rule call_peaks_seacr:
    input:
        treatment=OUTPUT_DIR + "/bedgraphs/{sample_id}.bedGraph",
    output:
        peaks=OUTPUT_DIR + "/peaks/seacr/{sample_id}.bed",
        seacr=temp(OUTPUT_DIR + "/peaks/seacr/{sample_id}_seacr.txt"),
        noM=temp(OUTPUT_DIR + "/bedgraphs/{sample_id}.nochrM.bedGraph"),
    wildcard_constraints:
        sample_id=r"(?!merged/).*",
    params:
        threshold=CONFIG.third_party_tools.seacr.threshold if CONFIG.third_party_tools.seacr else 0.01,
        norm=CONFIG.third_party_tools.seacr.normalization if CONFIG.third_party_tools.seacr else "non",
        stringency=CONFIG.third_party_tools.seacr.stringency if CONFIG.third_party_tools.seacr else "stringent",
        prefix=lambda wc, output: Path(output.peaks).parent
        / Path(output.peaks).name,
    threads: 1
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(
            initial_value=5, attempts=attempt, scale=SCALE_RESOURCES
        ),
        runtime=lambda wildcards, attempt: define_time_requested(
            initial_value=2, attempts=attempt, scale=SCALE_RESOURCES
        ),
    container:
        "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log:
        OUTPUT_DIR + "/logs/seacr/{sample_id}.log",
    benchmark:
        OUTPUT_DIR + "/.benchmark/seacr/{sample_id}.tsv"
    message:
        "Calling peaks with SEACR for sample {wildcards.sample_id}"
    shell:
        """
        awk '$1 != "chrM"' {input.treatment} > {output.noM}
        SEACR_1.3.sh {output.noM} {params.threshold} {params.norm} {params.stringency} {output.peaks} > {log} 2>&1 || touch {params.prefix}.{params.stringency}.bed
        mv {params.prefix}.{params.stringency}.bed {output.seacr}
        cut -f 1-3 {output.seacr} > {output.peaks}
        """
