"""Helper functions for CRISPR workflows."""

import json
import shlex


def get_cutadapt_adapter_args(wildcards, CONFIG, OUTPUT_DIR):
    """
    Build cutadapt adapter arguments from detected adapters.
    Falls back to config if detection fails or returns None.

    Args:
        wildcards: Snakemake wildcards object containing 'sample'.
        CONFIG: The configuration object.
        OUTPUT_DIR: The output directory path.

    Returns:
        str: The cutadapt command line arguments with adapter specifications.
    """
    adapter_file = OUTPUT_DIR + f"/resources/{wildcards.sample}_adapters.json"
    base_options = str(CONFIG.third_party_tools.cutadapt.command_line_arguments)

    try:
        with open(adapter_file, "r") as f:
            adapters = json.load(f)

        adapter_r1 = adapters.get("adapter_r1")
        adapter_r2 = adapters.get("adapter_r2")

        tokens = shlex.split(base_options)
        filtered_tokens = []
        skip_next = False
        for i, token in enumerate(tokens):
            if skip_next:
                skip_next = False
                continue
            if token in ["-g", "-G", "-a", "-A"]:
                skip_next = True
                continue
            if (
                token.startswith("-g")
                or token.startswith("-G")
                or token.startswith("-a")
                or token.startswith("-A")
            ):
                continue
            filtered_tokens.append(token)

        base_options = " ".join(filtered_tokens)

        adapter_args = ""
        if adapter_r1:
            adapter_args += f" -g '{adapter_r1}'"
        if adapter_r2:
            adapter_args += f" -G '{adapter_r2}'"

        return base_options + adapter_args

    except (FileNotFoundError, PermissionError, json.JSONDecodeError, KeyError):
        return base_options
