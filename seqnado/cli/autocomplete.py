"""Typer autocomplete functions for CLI arguments."""

from __future__ import annotations

from pathlib import Path
from typing import List, Optional


def _assay_names() -> List[str]:
    """Get list of available assay names for autocomplete."""
    from seqnado.inputs import Assay  # local import to keep CLI startup snappy

    return list(Assay.all_assay_clean_names())


def assay_autocomplete(_: str) -> List[str]:
    """Autocomplete assay names."""
    return _assay_names()


def fastq_autocomplete(incomplete: str) -> List[str]:
    """Provide FASTQ file suggestions without expensive globbing."""
    p = Path(incomplete or ".")
    base = p.parent if p.name else p
    pattern = p.name or "*"
    try:
        candidates = sorted(base.glob(pattern + "*.fastq.gz"))
    except Exception:
        candidates = []
    return [str(x) for x in candidates if x.is_file()][:40]


def _get_profile_name(fn: Path) -> Optional[str]:
    """Wrapper for get_profile_name from utils (for backward compatibility)."""
    from seqnado.utils import get_profile_name

    return get_profile_name(fn)


def _preset_profiles() -> dict:
    """Wrapper for get_preset_profiles from utils (for backward compatibility)."""
    from seqnado.utils import get_preset_profiles

    return get_preset_profiles()


def _profile_autocomplete() -> List[str]:
    """Return list of profile shortcode keys for CLI tab-completion."""
    return list(_preset_profiles().keys())


def _find_fastqs(hints: List[str]) -> List[Path]:
    """
    Search the provided hint directories for *.fastq.gz files.
    Skip hints that don't exist and return a sorted list.
    Searches both directly in the hint directory and in subdirectories.
    """
    seen: set[Path] = set()
    out: List[Path] = []
    for loc in hints:
        p = Path(loc)
        if not p.exists():
            continue
        # Search directly in the directory and one level deeper
        for pattern in ("*.fastq.gz", "*/*.fastq.gz"):
            for f in sorted(p.glob(pattern)):
                resolved = f.resolve()
                if resolved not in seen:
                    seen.add(resolved)
                    out.append(f)
    return out
