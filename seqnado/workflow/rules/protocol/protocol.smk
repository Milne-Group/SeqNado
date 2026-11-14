from seqnado._version import __version__ as SEQNADO_VERSION

assay_for_protocol = ASSAY.name

protocol_inputs = [
    str(p)
    for p in OUTPUT.files
    if "/geo_submission/" not in str(p) 
    and not str(p).endswith("/protocol.txt")
    and "genome_browser_plots" not in str(p)
    and "track_plots" not in str(p)
    and "heatmap" not in str(p)
    and "hub/" not in str(p)
]


rule protocol:
    input:
        protocol_inputs,
    output:
        OUTPUT_DIR + "/protocol.txt",
    params:
        assay=assay_for_protocol,
        seqnado_version=SEQNADO_VERSION,
        # Tools with their own containers; versions are queried from these.
        dedicated_containers={
            "bamnado": CONTAINERS.BAMNADO,
            "mccnado": CONTAINERS.MCCNADO,
            "quantnado": CONTAINERS.QUANTNADO,
        },
    container: CONTAINERS.SEQNADO_PIPELINE
    log: OUTPUT_DIR + "/logs/protocol.log",
    benchmark: OUTPUT_DIR + "/.benchmark/protocol.tsv",
    message: "Producing data processing protocol",
    script:
        "../../scripts/produce_data_processing_protocol.py"