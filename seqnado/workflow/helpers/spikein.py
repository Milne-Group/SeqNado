"""Helpers for classifying contigs as reference (endogenous) or spike-in (exogenous).

The regexes built here are restricted to a POSIX-ERE-compatible subset so the
same pattern can be used by ``awk`` inside Snakemake shell rules and by Python's
``re`` module in unit tests.
"""

from __future__ import annotations

import re
from collections.abc import Iterable, Sequence

# Contigs starting with "chr" (UCSC style, including scaffolds such as
# chr1_KI270706v1_random) or named like Ensembl/NCBI primary chromosomes
# (1, 2, ..., X, Y, Z, W, M, MT).
DEFAULT_REFERENCE_CONTIG_REGEX = r"^(chr|([0-9]+|[XYZW]|MT?)$)"

_ERE_SPECIAL = set(r".[]()*+?{}|^$\\")


def ere_escape(text: str) -> str:
    """Escape ``text`` so it matches literally in both POSIX ERE and Python ``re``."""
    return "".join(f"\\{c}" if c in _ERE_SPECIAL else c for c in text)


def _clean(prefixes: Iterable[str] | None) -> list[str]:
    if prefixes is None:
        return []
    if isinstance(prefixes, str):
        prefixes = [prefixes]
    return [p.strip() for p in prefixes if p and p.strip()]


def reference_contig_regex(prefixes: Sequence[str] | str | None = None) -> str:
    """Regex matching reference-genome contigs.

    If ``prefixes`` is given, a contig matches when it starts with any of them.
    Otherwise :data:`DEFAULT_REFERENCE_CONTIG_REGEX` is used.
    """
    cleaned = _clean(prefixes)
    if not cleaned:
        return DEFAULT_REFERENCE_CONTIG_REGEX
    return "^(" + "|".join(ere_escape(p) for p in cleaned) + ")"


def spikein_contig_regex(prefixes: Sequence[str] | str) -> str:
    """Regex matching spike-in contigs.

    A contig matches when it equals one of ``prefixes`` or starts with
    ``<prefix>_`` (the convention used by ``seqnado genomes build --spikein``).
    """
    cleaned = _clean(prefixes)
    if not cleaned:
        raise ValueError("At least one spike-in contig prefix is required.")
    return "^(" + "|".join(ere_escape(p) for p in cleaned) + ")(_|$)"


def spikein_regexes_from_config(spikein_config) -> tuple[str, str]:
    """Build ``(reference_regex, spikein_regex)`` from a ``SpikeInConfig``."""
    exo_prefixes = spikein_config.exogenous_contig_prefix or _clean(
        spikein_config.exogenous_genome
    )
    if not exo_prefixes:
        raise ValueError(
            "Spike-in splitting needs 'exogenous_genome' or 'exogenous_contig_prefix' "
            "set in the spikein section of the config."
        )
    return (
        reference_contig_regex(spikein_config.endogenous_contig_prefix),
        spikein_contig_regex(exo_prefixes),
    )


def classify_contigs(
    contigs: Iterable[str], reference_regex: str, spikein_regex: str
) -> tuple[list[str], list[str], list[str]]:
    """Split contig names into ``(reference, spikein, unassigned)``.

    Spike-in matches take precedence over reference matches, mirroring the
    shell implementation in ``split_bam_by_contig.sh``.
    """
    ref_re = re.compile(reference_regex)
    exo_re = re.compile(spikein_regex)
    ref, exo, unassigned = [], [], []
    for contig in contigs:
        if exo_re.search(contig):
            exo.append(contig)
        elif ref_re.search(contig):
            ref.append(contig)
        else:
            unassigned.append(contig)
    return ref, exo, unassigned
