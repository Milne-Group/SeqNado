from pathlib import Path
import re
import numpy as np
from seqnado import Assay
from seqnado.workflow.helpers.visualise import get_hub_input_files

# Only create hub if enabled AND bigwigs/bigbeds will be created
if CONFIG.assay_config.create_ucsc_hub and (OUTPUT.bigwig_files or OUTPUT.bigbed_files or OUTPUT.sentinel_files):
    
    rule generate_hub:
        input:
            data=lambda wc: get_hub_input_files(wc, OUTPUT, ASSAY, rules),
            report=OUTPUT_DIR + "/seqnado_report.html",
        output:
            hub=OUTPUT.ucsc_hub_files,
        params:
            assay=ASSAY,
            has_consensus_peaks=OUTPUT.has_consensus_peaks,
            genome=CONFIG.assay_config.ucsc_hub.genome,
            hub_name=CONFIG.assay_config.ucsc_hub.name,
            hub_email=CONFIG.assay_config.ucsc_hub.email,
            custom_genome=lambda wc: True if CONFIG.assay_config.ucsc_hub.two_bit else False,
            genome_twobit=CONFIG.assay_config.ucsc_hub.two_bit,
            genome_organism=CONFIG.assay_config.ucsc_hub.organism,
            genome_default_position=CONFIG.assay_config.ucsc_hub.default_position,
            color_by=CONFIG.assay_config.ucsc_hub.color_by,
            subgroup_by=CONFIG.assay_config.ucsc_hub.subgroup_by,
            supergroup_by=CONFIG.assay_config.ucsc_hub.supergroup_by,
            overlay_by=CONFIG.assay_config.ucsc_hub.overlay_by,
        containerized: False
        log: OUTPUT_DIR + f"/logs/{CONFIG.assay_config.ucsc_hub.name}.hub.log",
        benchmark: OUTPUT_DIR + f"/.benchmark/visualise/{CONFIG.assay_config.ucsc_hub.name}_hub.tsv",
        message: f"Generating UCSC Genome Browser hub: {CONFIG.assay_config.ucsc_hub.name}"
        script: "../../scripts/create_hub_with_tracknado.py"


    localrules:
        generate_hub