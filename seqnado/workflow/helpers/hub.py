"""Metadata extraction for UCSC hub generation (tracknado).

tracknado groups/colours tracks by columns of a metadata table. The built-in
``tracknado.extractors.from_seqnado_path`` assumes a legacy layout with an
``{assay}`` directory (``.../seqnado_output/{assay}/bigwigs/...``) that SeqNado
does not actually produce - real paths are ``seqnado_output/bigwigs/{method}/
{scale}/{sample}.bigWig`` with no assay directory, and spike-in/merged/
comparison variants add extra levels. Against real paths that extractor mis-
labels ``assay``/``method``/``norm`` and never emits ``antibody``/``strand``/
``viewpoint``, so the per-assay hub defaults crash with a missing-column
assertion inside ``TrackDesign``.

This module replaces it with an assay-aware extractor. The Snakemake rule
already knows the real assay (``snakemake.params.assay``), so we pass it in
rather than guessing it from the path. Every track of a given assay is returned
with the *same* set of keys (filling ``"NA"`` where a field does not apply), so
any column a hub config groups on is guaranteed to exist for all tracks.

Real path shapes handled (relative to ``seqnado_output``):

    bigwigs/{method}/{scale}/{sample}.bigWig                 # individual
    bigwigs/{method}/{scale}/{sample}_{plus|minus}.bigWig    # RNA, stranded
    bigwigs/{method}/spikein/{spikein}/{sample}.bigWig       # spike-in scaled
    bigwigs/{method}/merged/{scale}/{group}.bigWig           # merged replicates
    bigwigs/{method}/aggregated/{cond}.bigWig                # condition mean
    bigwigs/{method}/subtraction/{c1}_vs_{c2}.bigWig         # condition diff
    bigwigs/mcc/replicates/{sample}_{viewpoint}.bigWig       # MCC per viewpoint
    bigwigs/mcc/aggregated-using-mean/{group}_{viewpoint}.bigWig
    bigwigs/mcc/subtractions/{g1}_vs_{g2}_{viewpoint}.bigWig
    peaks/{method}/{sample}.bed  (and .../merged/{sample}.bed, -> .bb)
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Callable, Mapping, Sequence

import pandas as pd

from seqnado.inputs.core import Metadata

# Scaling techniques that appear as a directory level (DataScalingTechnique).
_SCALING_TOKENS = {"unscaled", "csaw", "cpm", "rpkm", "spikein"}
# Spike-in normalisation methods (SpikeInMethod) - the level after "spikein/".
_SPIKEIN_TOKENS = {"orlando", "with_input", "deseq2", "edger"}
# Directory levels that mark a condition-level comparison track.
_AGGREGATED_TOKENS = {"aggregated", "aggregated-using-mean"}
_SUBTRACTION_TOKENS = {"subtraction", "subtractions"}

_TRACK_SUFFIXES = (".bigWig", ".bigwig", ".bw", ".bb", ".bigBed", ".bed")

# Fields derived from the output *path* (present for every track).
_PATH_KEYS = (
    "assay",
    "file_type",
    "method",
    "norm",
    "spikein",
    "samplename",
    "strand",
    "viewpoint",
    "antibody",
    "merged",
    "track_kind",
)

# Fields pulled from the design table (metadata.csv), keyed by sample. These
# are the standard SeqNado metadata fields; ``load_design_metadata`` also
# retains user-defined columns so a hub can group by any column in the design.
# "assay" is dropped because the path-derived "assay" already carries it.
DESIGN_COLUMNS = tuple(name for name in Metadata.model_fields if name != "assay")

# Every extracted track carries the *same* set of keys regardless of assay, so a
# column referenced by a hub config is guaranteed present on all tracks
# (tracknado's composite/overlay grouping asserts this). Fields that do not
# apply to a given track are "NA".
METADATA_KEYS = _PATH_KEYS + DESIGN_COLUMNS


def load_design_metadata(path: str | Path | None) -> dict[str, dict[str, str]]:
    """Load per-sample design columns from a SeqNado metadata.csv.

    Returns a mapping ``sample_id -> {column: value}`` for every non-identifier
    column in the file. Missing/NA values become ``"NA"``. ``uid`` is also
    registered as an alias when present, which supports designs whose output
    sample names are based on that field. Path-derived metadata always takes
    precedence over a same-named design column.

    Returns an empty mapping if the file is absent or has no ``sample_id``
    column, so hub generation degrades gracefully.
    """
    if not path:
        return {}
    path = Path(path)
    if not path.exists():
        return {}

    df = pd.read_csv(path)
    if "sample_id" not in df.columns:
        return {}

    cols = [c for c in df.columns if c not in {"sample_id", "uid"}]
    mapping: dict[str, dict[str, str]] = {}
    for _, row in df.iterrows():
        sample_id = str(row["sample_id"])
        values = {
            c: ("NA" if pd.isna(row[c]) else str(row[c])) for c in cols
        }
        mapping[sample_id] = values
        if "uid" in df.columns and not pd.isna(row["uid"]):
            mapping[str(row["uid"])] = values
    return mapping


def _stem(path: Path) -> str:
    """Filename without its (possibly compound) track suffix."""
    name = path.name
    for suffix in _TRACK_SUFFIXES:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return path.stem


def extract_seqnado_track_metadata(
    path: Path,
    assay: str,
    design_metadata: Mapping[str, Mapping[str, str]] | None = None,
) -> dict[str, str]:
    """Extract hub grouping metadata from a SeqNado output path.

    Args:
        path: Path to a bigWig/bigBed/bed track produced by SeqNado.
        assay: The real assay clean name (e.g. ``"chip"``, ``"rna"``). Passed in
            from the workflow rather than parsed from the path.
        design_metadata: Optional mapping returned by :func:`load_design_metadata`.
            Its fields are made available for hub grouping and colouring.

    Returns:
        A dict with every key in :data:`METADATA_KEYS`; fields that do not apply
        to this track are ``"NA"``.
    """
    assay = assay.lower()
    parts = [p for p in path.parts]
    parts_lower = [p.lower() for p in parts]

    design_columns = tuple(
        sorted({column for values in (design_metadata or {}).values() for column in values})
    )
    meta: dict[str, str] = {k: "NA" for k in (*METADATA_KEYS, *design_columns)}
    meta["assay"] = assay

    # Locate the bigwigs/peaks anchor; the level after it is the method.
    anchor_idx = None
    for i, part in enumerate(parts_lower):
        if part in ("bigwigs", "peaks"):
            anchor_idx = i
            break

    if anchor_idx is None:
        # Unknown layout - fall back to a bare samplename so the hub still builds.
        meta["samplename"] = _stem(path)
        return _add_design_metadata(meta, design_metadata)

    meta["file_type"] = "signal" if parts_lower[anchor_idx] == "bigwigs" else "peaks"
    if anchor_idx + 1 < len(parts):
        meta["method"] = parts[anchor_idx + 1]

    # Directory levels between the method and the filename encode scale, spike-in
    # method, merged status and comparison kind (order varies, so scan them all).
    mids = parts_lower[anchor_idx + 2 : -1]
    meta["track_kind"] = "individual"
    seen_spikein = False
    for level in mids:
        if level == "spikein":
            meta["norm"] = "spikein"
            seen_spikein = True
        elif seen_spikein and meta["spikein"] == "NA" and level in _SPIKEIN_TOKENS:
            meta["spikein"] = level
        elif level in _SCALING_TOKENS:
            meta["norm"] = level
        elif level == "merged":
            meta["merged"] = "yes"
        elif level in _AGGREGATED_TOKENS:
            meta["track_kind"] = "aggregated"
        elif level in _SUBTRACTION_TOKENS:
            meta["track_kind"] = "subtraction"
        elif level == "replicates":
            meta["track_kind"] = "individual"

    # Parse the filename stem for assay-specific fields.
    stem = _stem(path)
    meta["samplename"] = stem
    is_individual = meta["track_kind"] == "individual"

    if assay == "rna":
        m = re.search(r"[._](plus|minus)$", stem)
        if m:
            meta["strand"] = m.group(1)
            meta["samplename"] = stem[: m.start()]
    elif assay == "mcc":
        # Every MCC track carries a trailing _{viewpoint}.
        if "_" in stem:
            meta["samplename"], meta["viewpoint"] = stem.rsplit("_", 1)
    elif assay in ("chip", "cat"):
        # Individual IP tracks are named {sample}_{antibody}; comparison tracks
        # (aggregated/subtraction) are condition-level and carry no antibody.
        if is_individual and "_" in stem:
            meta["samplename"], meta["antibody"] = stem.rsplit("_", 1)

    return _add_design_metadata(meta, design_metadata)


def _add_design_metadata(
    meta: dict[str, str],
    design_metadata: Mapping[str, Mapping[str, str]] | None,
) -> dict[str, str]:
    """Add design values for a track without overwriting path-derived fields."""
    if not design_metadata:
        return meta

    values = design_metadata.get(meta["samplename"])
    if not values:
        return meta

    for column, value in values.items():
        if column not in _PATH_KEYS:
            meta[column] = value
    return meta


def make_seqnado_extractor(
    assay: str,
    design_metadata: Mapping[str, Mapping[str, str]] | None = None,
) -> Callable[[Path], dict[str, str]]:
    """Return a tracknado metadata extractor bound to an assay and design."""

    def _extractor(path: Path) -> dict[str, str]:
        return extract_seqnado_track_metadata(Path(path), assay, design_metadata)

    return _extractor


def format_hub_labels(hub: Any, supergroup_by: Sequence[str] | None) -> None:
    """Give TrackNado-generated containers concise, human-readable UCSC labels.

    TrackNado uses its internal identifiers (for example,
    ``signal_treated_H3K27ac``) as both labels. UCSC limits ``shortLabel`` to
    17 characters, so those identifiers are truncated in the browser. Keep the
    stable identifiers unchanged, but replace their display labels after the
    hub is built and before it is staged.
    """
    if not supergroup_by or not getattr(hub.track_design, "super_tracks", None):
        return

    details = hub.track_design.details
    for supertrack_id, supertrack in hub.track_design.super_tracks.items():
        rows = details.loc[details["supertrack"] == supertrack_id]
        if rows.empty:
            continue

        row = rows.iloc[0]
        groups = [(column, str(row[column])) for column in supergroup_by]
        long_label = " · ".join(
            f"{column.replace('_', ' ').title()}: {value}"
            for column, value in groups
            if value != "NA"
        )
        short_values = [value for column, value in groups if column != "file_type"]
        short_label = _short_hub_label(" ".join(short_values) or long_label)
        _set_track_labels(supertrack, short_label, long_label)

        for child in supertrack.subtracks:
            child_type = "Signal" if child.tracktype == "bigWig" else "Peaks"
            _set_track_labels(
                child,
                _short_hub_label(f"{child_type} tracks"),
                f"{child_type} tracks — {long_label}",
            )


def _short_hub_label(value: str) -> str:
    """Return a UCSC-compliant label, preferring whole words when trimming."""
    limit = 17
    if len(value) <= limit:
        return value
    shortened = value[:limit].rsplit(" ", 1)[0]
    return shortened or value[:limit]


def _set_track_labels(track: Any, short_label: str, long_label: str) -> None:
    """Update trackhub's attributes and serialised keyword arguments together."""
    track.short_label = short_label
    track.long_label = long_label
    if hasattr(track, "add_params"):
        track.add_params(shortLabel=short_label, longLabel=long_label)
    else:
        track.kwargs["shortLabel"] = short_label
        track.kwargs["longLabel"] = long_label
