"""Helpers for QuantNado dataset rules."""

from typing import Any


def get_quantnado_pairedness_flag(input_files: Any, sample_id: str) -> str:
    """Return the QuantNado pairedness flag for a sample.

    QuantNado defaults BAM inputs to paired-end, so only single-end samples need
    an explicit flag.
    """
    is_paired_end = getattr(input_files, "is_paired_end", None)
    if not callable(is_paired_end):
        return ""

    try:
        return "" if is_paired_end(sample_id) else "--single-end"
    except (KeyError, TypeError, AttributeError):
        return ""
