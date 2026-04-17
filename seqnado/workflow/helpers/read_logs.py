"""Helpers for writing standardised read-count log TSV rows from Snakemake rules."""

from __future__ import annotations


READ_LOG_HEADER = "Step\\tSample\\tReads Before\\tReads After\\n"


def read_log_shared_path(output_dir: str, entity: str = "{sample}", subdir: str = "alignment_post_process") -> str:
    """Return the shared per-entity read log path."""
    return f"{output_dir}/qc/{subdir}/{entity}.tsv"


def emit_read_logs(
    step: str,
    sample_token: str,
    shared_log: str,
    read_log: str | None = None,
) -> str:
    """Return shell commands to append a standardised read-count TSV row."""
    row_format = f"{step}\\t%s\\t%s\\t%s\\n"
    commands = [
        f"mkdir -p $(dirname {shared_log})",
        f"if [ ! -f {shared_log} ]; then printf '{READ_LOG_HEADER}' > {shared_log}; fi",
    ]
    if read_log is not None:
        commands.extend(
            [
                f"mkdir -p $(dirname {read_log})",
                f"printf '{READ_LOG_HEADER}{row_format}' '{sample_token}' \"$before\" \"$after\" > {read_log}",
            ]
        )
    commands.append(
        f"printf '{row_format}' '{sample_token}' \"$before\" \"$after\" >> {shared_log}"
    )
    return " &&\n    ".join(commands)
