"""Post-pipeline analysis API for SeqNado output directories.

Usage:

    from seqnado.analysis import SeqNadoProject

    proj = SeqNadoProject("seqnado_output")
    # or with explicit config
    proj = SeqNadoProject("seqnado_output", config="chip_config.yml")

    proj.samples          # list[str]
    proj.conditions       # list[str]
    proj.antibodies       # list[str]  (ChIP / CUT&Tag)
    proj.design           # pd.DataFrame  (read/write — see below)

    # File paths (list[Path], filtered to existing files)
    proj.bigwigs(condition="treated", method="deeptools", scale="unscaled")
    proj.bigwigs(antibody="H3K27ac", merged=True)
    proj.peaks(method="macs2", merged=False)
    proj.bams(antibody="H3K27ac")
    proj.track_plots(method="deeptools", scale="unscaled")
    proj.contacts()
    proj.logs(rule="bowtie2", sample="sample1")

    # Load data as DataFrames
    proj.load_counts()                 # featureCounts / Salmon matrix
    proj.load_peaks()                  # all BED peaks concatenated
    proj.load_alignment_stats()        # per-sample alignment step counts
    proj.load_frip()                   # FRiP scores per sample / peak method
    proj.load_library_complexity()     # Picard duplicate metrics
    proj.load_benchmarks()             # Snakemake rule runtimes
    proj.load_normalisation_factors()  # spike-in / CSAW factors

    # Convenience
    proj.protocol()                    # Path to protocol.txt
    proj.bai_for(bam_path)             # BAI index for a BAM
    proj.sample_files("sample1")       # dict of all files for one sample
    proj.condition_pairs()             # [("ctrl", "treated"), ...]

    # Filtered views (chainable)
    proj.filter(condition="treated").bigwigs()
    proj.filter(antibody="H3K27ac").peaks(merged=True)

    # DataFrame indexes with full metadata
    proj.bigwig_dataframe()
    proj.peak_dataframe()

    # Updating metadata
    # ------------------
    # ``proj.design`` is writable.  Assigning a new DataFrame clears all cached
    # metadata (conditions, antibodies, groups, _sample_meta) so that every
    # subsequent access re-derives from the new data.  Use this to add or
    # correct columns that are absent from the original design file.
    #
    # Example — add conditions by stripping the replicate suffix from sample IDs:
    #
    #     df = proj.design.copy()
    #     df["condition"] = df["sample_id"].str.replace(r"-\\d+$", "", regex=True)
    #     proj.design = df
    #     proj.conditions   # now populated
"""

from __future__ import annotations

from enum import Enum
from functools import cached_property
from pathlib import Path
from typing import TYPE_CHECKING, Union

import pandas as pd

from seqnado.core import DataScalingTechnique, PeakCallingMethod, PileupMethod, SpikeInMethod

if TYPE_CHECKING:
    from seqnado.config import SeqnadoConfig


# ---------------------------------------------------------------------------
# Type aliases
# ---------------------------------------------------------------------------

# Enum-or-string convenience type for filter parameters
_MethodArg = Union[str, PileupMethod]
_ScaleArg = Union[str, DataScalingTechnique]
_PeakMethodArg = Union[str, PeakCallingMethod]
_SpikeInArg = Union[str, SpikeInMethod]


def _ev(v: Enum | str | None) -> str | None:
    """Extract string value from an enum, or return the string unchanged."""
    if v is None:
        return None
    return v.value if isinstance(v, Enum) else v


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _load_yaml(path: Path) -> dict:
    try:
        import yaml
    except ImportError:
        raise ImportError("PyYAML is required: pip install pyyaml")
    with open(path) as fh:
        return yaml.safe_load(fh)


def _parse_bigwig_path(p: Path, bigwig_dir: Path) -> dict | None:
    """Parse a bigWig path into metadata fields.

    Expected structures (relative to bigwig_dir):
      {method}/{scale}/{sample}.bigWig
      {method}/merged/{scale}/{sample}.bigWig
      {method}/spikein/{spikein_method}/{sample}.bigWig
      {method}/merged/spikein/{spikein_method}/{sample}.bigWig
      {method}/aggregated/{condition}.bigWig
      {method}/subtraction/{label}.bigWig
      {method}/{sample}_{genome}.bigWig           (methylation)
    """
    try:
        rel = p.relative_to(bigwig_dir)
    except ValueError:
        return None

    dir_parts = list(rel.parts[:-1])
    if not dir_parts:
        return None

    method = dir_parts[0]
    merged = "merged" in dir_parts
    spikein_method: str | None = None

    if "spikein" in dir_parts:
        idx = dir_parts.index("spikein")
        scale = "spikein"
        spikein_method = dir_parts[idx + 1] if idx + 1 < len(dir_parts) else None
    elif merged:
        merged_idx = dir_parts.index("merged")
        scale = dir_parts[merged_idx + 1] if merged_idx + 1 < len(dir_parts) else "unscaled"
    elif len(dir_parts) >= 2:
        scale = dir_parts[1]
    else:
        scale = "unscaled"

    stem = p.stem
    strand: str | None = None
    if stem.endswith("_plus"):
        strand = "plus"
        sample = stem[:-5]
    elif stem.endswith("_minus"):
        strand = "minus"
        sample = stem[:-6]
    else:
        sample = stem

    return {
        "path": p,
        "sample": sample,
        "method": method,
        "scale": scale,
        "spikein_method": spikein_method,
        "merged": merged,
        "strand": strand,
    }


def _parse_peak_path(p: Path, peaks_dir: Path) -> dict | None:
    try:
        rel = p.relative_to(peaks_dir)
    except ValueError:
        return None

    dir_parts = list(rel.parts[:-1])
    if not dir_parts:
        return None

    method = dir_parts[0]
    merged = "merged" in dir_parts
    return {"path": p, "sample": p.stem, "method": method, "merged": merged}


# ---------------------------------------------------------------------------
# Core class
# ---------------------------------------------------------------------------

class SeqNadoProject:
    """Analysis interface for a completed SeqNado output directory.

    Parameters
    ----------
    output_dir:
        Path to the SeqNado output directory (e.g. ``seqnado_output``).
    config:
        Path to the SeqNado YAML config file, or a pre-loaded
        ``SeqnadoConfig`` instance.  When omitted, the project directory
        is scanned for ``*_config.yml`` files automatically.
    design:
        Path to the design TSV / CSV, or a ``pd.DataFrame``.  When omitted
        the path is read from the config (``config.metadata``) or
        auto-discovered next to the output directory.
    """

    def __init__(
        self,
        output_dir: str | Path,
        config: str | Path | "SeqnadoConfig" | None = None,
        design: str | Path | pd.DataFrame | None = None,
    ):
        self._output_dir = Path(output_dir)
        self._config: "SeqnadoConfig | None" = self._load_config(config)
        self._design_df: pd.DataFrame | None = self._load_design(design)

    # ------------------------------------------------------------------ #
    #  Internal loaders                                                    #
    # ------------------------------------------------------------------ #

    def _load_config(self, config) -> "SeqnadoConfig | None":
        from seqnado.config import SeqnadoConfig as _SeqnadoConfig

        if isinstance(config, _SeqnadoConfig):
            self._config_path: Path | None = None
            return config

        path: Path | None = None
        if config is not None:
            path = Path(config)
        else:
            # CLI writes configs as config_{assay}.yaml; also accept legacy *_config.yml
            patterns = ["config_*.yaml", "config_*.yml", "*_config.yml"]
            for search_dir in [
                self._output_dir.parent,
                self._output_dir.parent.parent,
                Path("."),
            ]:
                for pattern in patterns:
                    candidates = sorted(search_dir.glob(pattern))
                    if candidates:
                        path = candidates[0]
                        break
                if path is not None:
                    break

        self._config_path = path if (path is not None and path.exists()) else None

        if self._config_path is None:
            return None

        try:
            data = _load_yaml(self._config_path)
            return _SeqnadoConfig.model_validate(
                data, context={"skip_path_validation": True}
            )
        except Exception:
            return None

    def _load_design(self, design) -> pd.DataFrame | None:
        if isinstance(design, pd.DataFrame):
            return design

        path: Path | None = None
        if design is not None:
            path = Path(design)
        elif self._config is not None and self._config.metadata is not None:
            # Primary source: metadata path recorded in the config YAML.
            # Resolve relative paths against the config file's directory.
            config_meta = Path(self._config.metadata)
            if not config_meta.is_absolute() and self._config_path is not None:
                config_meta = self._config_path.parent / config_meta
            if config_meta.exists():
                path = config_meta
        if path is None:
            # Assay name used to find multiomics-style metadata_{assay}.csv
            assay_name = (
                self._config.assay.clean_name
                if self._config is not None
                else self._output_dir.name
            )
            for candidate in [
                self._output_dir.parent / "design.tsv",
                self._output_dir.parent / "design.csv",
                self._output_dir / "design.tsv",
                self._output_dir / "design.csv",
                # project root when output_dir is a multiomics assay sub-dir
                self._output_dir.parent.parent / "design.tsv",
                self._output_dir.parent.parent / "design.csv",
                self._output_dir.parent.parent / f"metadata_{assay_name}.csv",
                self._output_dir.parent / f"metadata_{assay_name}.csv",
            ]:
                if candidate.exists():
                    path = candidate
                    break

        if path is None or not path.exists():
            return None

        sep = "," if path.suffix == ".csv" else "\t"
        df = pd.read_csv(path, sep=sep)

        # Normalise non-standard column aliases to the canonical names defined in
        # seqnado.inputs.validation.DesignDataFrame / seqnado.inputs.core.Metadata
        _ALIASES = {"sample_name": "sample_id", "scale_group": "scaling_group"}
        df = df.rename(columns={k: v for k, v in _ALIASES.items() if k in df.columns})
        return df

    # ------------------------------------------------------------------ #
    #  Indexes (lazily built, cached)                                      #
    # ------------------------------------------------------------------ #

    @cached_property
    def _bw_index(self) -> pd.DataFrame:
        bigwig_dir = self._output_dir / "bigwigs"
        if not bigwig_dir.exists():
            return pd.DataFrame(
                columns=["path", "sample", "method", "scale", "spikein_method", "merged", "strand"]
            )
        records = []
        for p in bigwig_dir.rglob("*.bigWig"):
            rec = _parse_bigwig_path(p, bigwig_dir)
            if rec:
                records.append(rec)
        return pd.DataFrame(records) if records else pd.DataFrame(
            columns=["path", "sample", "method", "scale", "spikein_method", "merged", "strand"]
        )

    @cached_property
    def _peak_index(self) -> pd.DataFrame:
        peaks_dir = self._output_dir / "peaks"
        if not peaks_dir.exists():
            return pd.DataFrame(columns=["path", "sample", "method", "merged"])
        records = []
        for p in peaks_dir.rglob("*.bed"):
            rec = _parse_peak_path(p, peaks_dir)
            if rec:
                records.append(rec)
        return pd.DataFrame(records) if records else pd.DataFrame(
            columns=["path", "sample", "method", "merged"]
        )

    # ------------------------------------------------------------------ #
    #  Sample ↔ metadata mappings                                          #
    # ------------------------------------------------------------------ #

    @cached_property
    def _sample_meta(self) -> pd.DataFrame:
        """DataFrame indexed by sample name with condition / antibody columns."""
        from seqnado.inputs.core import Metadata

        if self._design_df is None:
            return pd.DataFrame(columns=["condition", "antibody", "group"])

        df = self._design_df.copy()

        # ip is the canonical name in DesignDataFrame; rename to antibody for internal use
        if "ip" in df.columns and "antibody" not in df.columns:
            df = df.rename(columns={"ip": "antibody"})

        if "sample_id" not in df.columns:
            return pd.DataFrame(columns=["condition", "antibody", "group"])

        _meta_fields = set(Metadata.model_fields)
        rows = []
        for _, row in df.iterrows():
            sid = str(row["sample_id"])

            # Use Metadata model for proper coercion of None-like values
            meta = Metadata.model_validate(
                {k: row[k] for k in _meta_fields if k in df.columns}
            )

            ab_raw = row.get("antibody") if "antibody" in df.columns else None
            ab = None if pd.isna(ab_raw) else ab_raw

            rows.append({"sample": sid, "condition": meta.condition, "antibody": ab, "group": meta.scaling_group})

            if ab:
                rows.append({"sample": f"{sid}_{ab}", "condition": meta.condition, "antibody": ab, "group": meta.scaling_group})

        return pd.DataFrame(rows).drop_duplicates("sample").set_index("sample")

    def _resolve_samples(
        self,
        sample: str | list[str] | None,
        condition: str | list[str] | None,
        antibody: str | list[str] | None,
        group: str | list[str] | None,
    ) -> set[str] | None:
        """Return sample name set matching all provided filters, or None for no filter."""
        targets: set[str] | None = None

        if sample is not None:
            sl = [sample] if isinstance(sample, str) else list(sample)
            targets = set(sl)

        def _apply(col: str, values: str | list[str]) -> set[str]:
            vl = [values] if isinstance(values, str) else list(values)
            meta = self._sample_meta
            if col not in meta.columns:
                return set()
            matched = meta[meta[col].isin(vl)].index.tolist()
            return set(matched)

        for col, val in [("condition", condition), ("antibody", antibody), ("group", group)]:
            if val is not None:
                s = _apply(col, val)
                targets = s if targets is None else targets & s

        return targets

    # ------------------------------------------------------------------ #
    #  Public properties                                                   #
    # ------------------------------------------------------------------ #

    @property
    def output_dir(self) -> Path:
        return self._output_dir

    @property
    def config(self) -> "SeqnadoConfig | None":
        return self._config

    @property
    def design(self) -> pd.DataFrame | None:
        return self._design_df

    @design.setter
    def design(self, df: pd.DataFrame) -> None:
        """Replace the design DataFrame and invalidate all derived metadata.

        Use this to add or correct columns (e.g. ``condition``) that are
        absent from the original file::

            df = proj.design.copy()
            df["condition"] = df["sample_id"].str.replace(r"-\\d+$", "", regex=True)
            proj.design = df
        """
        self._design_df = df
        self.reload()

    @cached_property
    def samples(self) -> list[str]:
        """All sample names inferred from aligned BAMs."""
        bam_dir = self._output_dir / "aligned"
        if bam_dir.exists():
            return sorted(
                p.stem for p in bam_dir.glob("*.bam") if not p.name.endswith(".bai")
            )
        if self._design_df is not None and "sample_id" in self._design_df.columns:
            return sorted(self._design_df["sample_id"].tolist())
        return []

    @cached_property
    def conditions(self) -> list[str]:
        if self._design_df is None or "condition" not in self._design_df.columns:
            return []
        return sorted(self._design_df["condition"].dropna().unique().tolist())

    @cached_property
    def antibodies(self) -> list[str]:
        for col in ("ip", "antibody"):
            if self._design_df is not None and col in self._design_df.columns:
                return sorted(self._design_df[col].dropna().unique().tolist())
        return []

    @cached_property
    def groups(self) -> list[str]:
        if self._design_df is None or "scaling_group" not in self._design_df.columns:
            return []
        return sorted(self._design_df["scaling_group"].dropna().unique().tolist())

    @cached_property
    def assay(self) -> str | None:
        if self._config is not None:
            return self._config.assay.clean_name
        for p in self._output_dir.glob("seqnado_report_*.html"):
            return p.stem.removeprefix("seqnado_report_")
        return None

    # ------------------------------------------------------------------ #
    #  Discovery helpers                                                   #
    # ------------------------------------------------------------------ #

    @cached_property
    def pileup_methods(self) -> list[str]:
        """BigWig pileup methods present in the output (e.g. ``["deeptools", "bamnado"]``)."""
        if self._bw_index.empty:
            return []
        return sorted(self._bw_index["method"].dropna().unique().tolist())

    @cached_property
    def scales(self) -> list[str]:
        """Scaling methods present across all bigWigs (e.g. ``["unscaled", "csaw"]``)."""
        if self._bw_index.empty:
            return []
        return sorted(self._bw_index["scale"].dropna().unique().tolist())

    @cached_property
    def spikein_methods(self) -> list[str]:
        """Spike-in normalisation methods present (e.g. ``["orlando"]``)."""
        if self._bw_index.empty:
            return []
        vals = self._bw_index["spikein_method"].dropna().unique().tolist()
        return sorted(vals)

    @cached_property
    def peak_methods(self) -> list[str]:
        """Peak-calling methods present (e.g. ``["macs2", "homer"]``)."""
        if self._peak_index.empty:
            return []
        return sorted(self._peak_index["method"].dropna().unique().tolist())

    @cached_property
    def consensus_groups(self) -> list[str]:
        """Merged / consensus group names found in bigWig or peak outputs."""
        groups: set[str] = set()
        if not self._bw_index.empty:
            merged_bws = self._bw_index[self._bw_index["merged"]]
            groups.update(merged_bws["sample"].dropna().unique().tolist())
        if not self._peak_index.empty:
            merged_peaks = self._peak_index[self._peak_index["merged"]]
            groups.update(merged_peaks["sample"].dropna().unique().tolist())
        return sorted(groups)

    @cached_property
    def log_rules(self) -> list[str]:
        """Rule names that have log files present under ``logs/``."""
        log_dir = self._output_dir / "logs"
        if not log_dir.exists():
            return []
        return sorted(p.name for p in log_dir.iterdir() if p.is_dir())

    @cached_property
    def benchmark_rules(self) -> list[str]:
        """Rule names that have benchmark files present under ``.benchmark/``."""
        bench_dir = self._output_dir / ".benchmark"
        if not bench_dir.exists():
            return []
        return sorted(p.name for p in bench_dir.iterdir() if p.is_dir())

    def samples_for(
        self,
        condition: str | list[str] | None = None,
        antibody: str | list[str] | None = None,
        group: str | list[str] | None = None,
    ) -> list[str]:
        """Return sample names matching the given metadata filters.

        Equivalent to ``proj.filter(...).samples`` but without building
        a filtered view object.

        Examples
        --------
        ::

            proj.samples_for(condition="treated")
            proj.samples_for(antibody="H3K27ac")
            proj.samples_for(condition="treated", antibody="H3K27ac")
        """
        targets = self._resolve_samples(None, condition, antibody, group)
        if targets is None:
            return self.samples
        return sorted(targets & set(self.samples))

    def metadata_for(self, sample: str) -> dict:
        """Return metadata dict for a single sample (condition, antibody, group).

        Returns an empty dict if the sample is not in the design.
        """
        meta = self._sample_meta
        if sample not in meta.index:
            return {}
        row = meta.loc[sample]
        # Ensure scalar row (unique index guaranteed by _sample_meta build)
        if isinstance(row, pd.DataFrame):
            row = row.iloc[0]
        return {k: v for k, v in row.items() if not pd.isna(v)}

    def what_exists(self) -> pd.DataFrame:
        """Per-sample availability table.

        Returns a DataFrame with one row per sample (inferred from BAMs and
        bigWig index) and boolean columns for each major file type.
        Includes both individual and merged samples.
        """
        aligned_dir = self._output_dir / "aligned"
        bam_set: set[str] = set()
        if aligned_dir.exists():
            for p in aligned_dir.glob("*.bam"):
                if not p.name.endswith(".bai"):
                    bam_set.add(p.stem)
            merged_dir = aligned_dir / "merged"
            if merged_dir.exists():
                for p in merged_dir.glob("*.bam"):
                    if not p.name.endswith(".bai"):
                        bam_set.add(p.stem)

        bw_samples = set(self._bw_index["sample"].tolist()) if not self._bw_index.empty else set()
        peak_samples = set(self._peak_index["sample"].tolist()) if not self._peak_index.empty else set()

        meth_dir = self._output_dir / "methylation" / "methyldackel"
        meth_samples: set[str] = set()
        if meth_dir.exists():
            for p in meth_dir.glob("*.bedGraph"):
                meth_samples.add(p.stem.split("_")[0])

        vcf_dir = self._output_dir / "variant"
        vcf_samples: set[str] = set()
        if vcf_dir.exists():
            vcf_samples = {
                p.name.replace(".vcf.gz", "").replace(".anno", "")
                for p in vcf_dir.glob("*.vcf.gz")
            }

        all_samples = sorted(bam_set | bw_samples | peak_samples)
        rows = []
        for s in all_samples:
            meta = self.metadata_for(s)
            rows.append({
                "sample": s,
                "condition": meta.get("condition"),
                "antibody": meta.get("antibody"),
                "group": meta.get("group"),
                "bam": s in bam_set,
                "bigwig": s in bw_samples,
                "peaks": s in peak_samples,
                "methylation": s in meth_samples,
                "vcf": s in vcf_samples,
            })
        return pd.DataFrame(rows)

    # ------------------------------------------------------------------ #
    #  BigWig access                                                        #
    # ------------------------------------------------------------------ #

    def bigwig_dataframe(self) -> pd.DataFrame:
        """Full bigWig index as a DataFrame (path, sample, method, scale, …)."""
        return self._bw_index.copy()

    def bigwigs(
        self,
        *,
        sample: str | list[str] | None = None,
        condition: str | list[str] | None = None,
        antibody: str | list[str] | None = None,
        group: str | list[str] | None = None,
        method: _MethodArg | None = None,
        scale: _ScaleArg | None = None,
        spikein_method: _SpikeInArg | None = None,
        merged: bool | None = None,
        strand: str | None = None,
        only_existing: bool = True,
    ) -> list[Path]:
        """Return bigWig paths matching the given filters.

        Parameters
        ----------
        sample:      Sample name(s).
        condition:   Condition label(s) from the design file.
        antibody:    IP antibody / target name(s) (ChIP / CUT&Tag).
        group:       Scaling group name(s).
        method:      :class:`~seqnado.core.PileupMethod` enum or string,
                     e.g. ``PileupMethod.DEEPTOOLS`` or ``"deeptools"``.
        scale:       :class:`~seqnado.core.DataScalingTechnique` enum or string,
                     e.g. ``DataScalingTechnique.UNSCALED`` or ``"unscaled"``.
        spikein_method: :class:`~seqnado.core.SpikeInMethod` enum or string,
                     e.g. ``SpikeInMethod.ORLANDO`` or ``"orlando"``.
        merged:      ``True`` = consensus merged tracks only.
        strand:      ``"plus"`` or ``"minus"`` (RNA-seq).
        only_existing: Skip paths that don't exist on disk (default ``True``).
        """
        df = self._bw_index
        if df.empty:
            return []

        targets = self._resolve_samples(sample, condition, antibody, group)
        if targets is not None:
            df = df[df["sample"].isin(targets)]

        if method is not None:
            df = df[df["method"] == _ev(method)]
        if scale is not None:
            df = df[df["scale"] == _ev(scale)]
        if spikein_method is not None:
            df = df[df["spikein_method"] == _ev(spikein_method)]
        if merged is not None:
            df = df[df["merged"] == merged]
        if strand is not None:
            df = df[df["strand"] == strand]

        paths: list[Path] = df["path"].tolist()
        if only_existing:
            paths = [p for p in paths if p.exists()]
        return sorted(paths)

    # ------------------------------------------------------------------ #
    #  Peak access                                                         #
    # ------------------------------------------------------------------ #

    def peak_dataframe(self) -> pd.DataFrame:
        """Full peak index as a DataFrame (path, sample, method, merged)."""
        return self._peak_index.copy()

    def peaks(
        self,
        *,
        sample: str | list[str] | None = None,
        condition: str | list[str] | None = None,
        antibody: str | list[str] | None = None,
        group: str | list[str] | None = None,
        method: _PeakMethodArg | None = None,
        merged: bool | None = None,
        only_existing: bool = True,
    ) -> list[Path]:
        """Return peak BED paths matching the given filters.

        ``method`` accepts a :class:`~seqnado.core.PeakCallingMethod` enum or
        string, e.g. ``PeakCallingMethod.MACS2`` or ``"macs2"``.
        """
        df = self._peak_index
        if df.empty:
            return []

        targets = self._resolve_samples(sample, condition, antibody, group)
        if targets is not None:
            df = df[df["sample"].isin(targets)]
        if method is not None:
            df = df[df["method"] == _ev(method)]
        if merged is not None:
            df = df[df["merged"] == merged]

        paths: list[Path] = df["path"].tolist()
        if only_existing:
            paths = [p for p in paths if p.exists()]
        return sorted(paths)

    # ------------------------------------------------------------------ #
    #  BAM access                                                          #
    # ------------------------------------------------------------------ #

    def bams(
        self,
        *,
        sample: str | list[str] | None = None,
        condition: str | list[str] | None = None,
        antibody: str | list[str] | None = None,
        group: str | list[str] | None = None,
        merged: bool = False,
        only_existing: bool = True,
    ) -> list[Path]:
        """Return aligned BAM paths."""
        bam_dir = self._output_dir / "aligned"
        if merged:
            bam_dir = bam_dir / "merged"

        if not bam_dir.exists():
            return []

        all_bams: dict[str, Path] = {
            p.stem: p
            for p in bam_dir.glob("*.bam")
            if not p.name.endswith(".bai")
        }

        targets = self._resolve_samples(sample, condition, antibody, group)
        if targets is not None:
            paths = [all_bams[s] for s in sorted(targets) if s in all_bams]
        else:
            paths = sorted(all_bams.values())

        if only_existing:
            paths = [p for p in paths if p.exists()]
        return paths

    # ------------------------------------------------------------------ #
    #  Counts / quantification                                             #
    # ------------------------------------------------------------------ #

    def counts(self, method: str = "feature_counts") -> Path | None:
        """Return path to the counts matrix for the given quantification method.

        Methods: ``"feature_counts"`` (default), ``"salmon"``.
        Returns ``None`` when the file does not exist.
        """
        paths = {
            "feature_counts": self._output_dir / "readcounts" / "feature_counts" / "read_counts.tsv",
            "salmon": self._output_dir / "readcounts" / "salmon" / "salmon_counts.csv",
        }
        if method not in paths:
            raise ValueError(f"Unknown method '{method}'. Choose from: {list(paths)}")
        p = paths[method]
        return p if p.exists() else None

    def load_counts(self, method: str = "feature_counts", **read_kwargs) -> pd.DataFrame:
        """Load the counts matrix as a ``pd.DataFrame``.

        Extra keyword args are forwarded to ``pd.read_csv``.
        """
        p = self.counts(method)
        if p is None:
            raise FileNotFoundError(f"No counts file found for method '{method}'")
        sep = "\t" if p.suffix == ".tsv" else ","
        read_kwargs.setdefault("sep", sep)
        read_kwargs.setdefault("index_col", 0)
        return pd.read_csv(p, **read_kwargs)

    # ------------------------------------------------------------------ #
    #  QC / reports                                                        #
    # ------------------------------------------------------------------ #

    def qc(
        self,
        *,
        sample: str | list[str] | None = None,
        tool: str | None = None,
        only_existing: bool = True,
    ) -> list[Path]:
        """Return QC report HTML files.

        Parameters
        ----------
        sample: Filter to a specific sample name.
        tool:   Filter by tool directory name, e.g. ``"fastqc_raw"``,
                ``"qualimap_bamqc"``, ``"fastq_screen"``.
        """
        qc_dir = self._output_dir / "qc"
        if not qc_dir.exists():
            return []

        search_dir = qc_dir / tool if tool else qc_dir
        paths = sorted(search_dir.rglob("*.html"))

        if sample is not None:
            sl = [sample] if isinstance(sample, str) else list(sample)
            paths = [p for p in paths if any(s in p.name for s in sl)]

        if only_existing:
            paths = [p for p in paths if p.exists()]
        return paths

    def report(self) -> Path | None:
        """Path to the MultiQC summary report HTML."""
        for p in sorted(self._output_dir.glob("seqnado_report_*.html")):
            return p
        return None

    # ------------------------------------------------------------------ #
    #  Assay-specific                                                      #
    # ------------------------------------------------------------------ #

    def methylation(
        self,
        *,
        sample: str | list[str] | None = None,
        condition: str | list[str] | None = None,
        genome: str | None = None,
        inverted: bool | None = None,
        only_existing: bool = True,
    ) -> list[Path]:
        """Return MethylDackel bedGraph files.

        Parameters
        ----------
        genome:   Filter to a specific reference genome name, e.g. ``"hg38"``.
        inverted: ``True`` = TAPS-inverted files only; ``False`` = regular CpG files only.
        """
        meth_dir = self._output_dir / "methylation" / "methyldackel"
        if not meth_dir.exists():
            return []

        paths = list(meth_dir.glob("*.bedGraph"))

        targets = self._resolve_samples(sample, condition, None, None)
        if targets is not None:
            paths = [p for p in paths if any(p.name.startswith(f"{t}_") for t in targets)]

        if genome is not None:
            paths = [p for p in paths if f"_{genome}_" in p.name or p.name.endswith(f"_{genome}.bedGraph")]

        if inverted is True:
            paths = [p for p in paths if p.name.endswith("_inverted.bedGraph")]
        elif inverted is False:
            paths = [p for p in paths if not p.name.endswith("_inverted.bedGraph")]

        if only_existing:
            paths = [p for p in paths if p.exists()]
        return sorted(paths)

    def vcf(
        self,
        *,
        sample: str | list[str] | None = None,
        annotated: bool = False,
        only_existing: bool = True,
    ) -> list[Path]:
        """Return variant VCF paths."""
        variant_dir = self._output_dir / "variant"
        if not variant_dir.exists():
            return []

        if annotated:
            vcfs: dict[str, Path] = {
                p.name.replace(".anno.vcf.gz", ""): p
                for p in variant_dir.glob("*.anno.vcf.gz")
            }
        else:
            vcfs = {
                p.name.replace(".vcf.gz", ""): p
                for p in variant_dir.glob("*.vcf.gz")
                if ".anno." not in p.name
            }

        if sample is not None:
            sl = [sample] if isinstance(sample, str) else list(sample)
            vcfs = {k: v for k, v in vcfs.items() if k in sl}

        paths = sorted(vcfs.values())
        if only_existing:
            paths = [p for p in paths if p.exists()]
        return paths

    def heatmaps(
        self,
        *,
        method: _MethodArg | None = None,
        scale: _ScaleArg | None = None,
        merged: bool | None = None,
        plot_type: str = "heatmap",
        only_existing: bool = True,
    ) -> list[Path]:
        """Return heatmap PDF paths.

        Parameters
        ----------
        method:    :class:`~seqnado.core.PileupMethod` enum or string.
        scale:     :class:`~seqnado.core.DataScalingTechnique` enum or string.
        plot_type: ``"heatmap"`` (default) or ``"metaplot"``.
        """
        heatmap_dir = self._output_dir / "heatmap"
        if not heatmap_dir.exists():
            return []

        paths = list(heatmap_dir.rglob(f"{plot_type}.pdf"))

        method_s = _ev(method)
        scale_s = _ev(scale)
        if method_s is not None:
            paths = [p for p in paths if f"/{method_s}/" in p.as_posix()]
        if scale_s is not None:
            paths = [p for p in paths if f"/{scale_s}/" in p.as_posix() or p.parent.name == scale_s]
        if merged is True:
            paths = [p for p in paths if "/merged/" in p.as_posix()]
        elif merged is False:
            paths = [p for p in paths if "/merged/" not in p.as_posix()]

        if only_existing:
            paths = [p for p in paths if p.exists()]
        return sorted(paths)

    def normalisation_factors(
        self, method: _SpikeInArg | None = None
    ) -> list[Path]:
        """Return normalisation factor TSV files.

        ``method`` accepts a :class:`~seqnado.core.SpikeInMethod` enum or string.
        """
        resources_dir = self._output_dir / "resources"
        if not resources_dir.exists():
            return []
        if method is not None:
            method_s = _ev(method) or ""
            p = resources_dir / method_s / "normalisation_factors.tsv"
            return [p] if p.exists() else []
        return sorted(resources_dir.rglob("normalisation_factors.tsv"))

    def load_spikein_stats(self) -> pd.DataFrame:
        """Load per-sample spike-in alignment counts as a ``pd.DataFrame``.

        Reads all ``aligned/spikein/{sample}_stats.tsv`` files produced by
        the ``bam_split`` rule.  Columns: ``sample``, ``reference_reads``,
        ``spikein_reads``.
        """
        spikein_dir = self._output_dir / "aligned" / "spikein"
        if not spikein_dir.exists():
            raise FileNotFoundError(f"Spike-in directory not found: {spikein_dir}")

        frames = []
        for p in sorted(spikein_dir.glob("*_stats.tsv")):
            frames.append(pd.read_csv(p, sep="\t"))

        if not frames:
            raise FileNotFoundError("No spike-in stats files found.")
        return pd.concat(frames, ignore_index=True)

    def load_normalisation_factors(self, method: _SpikeInArg | None = None) -> pd.DataFrame:
        """Load normalisation factors as a ``pd.DataFrame``."""
        paths = self.normalisation_factors(method)
        if not paths:
            raise FileNotFoundError("No normalisation factor files found.")
        frames = []
        for p in paths:
            df = pd.read_csv(p, sep="\t")
            df["method"] = p.parent.name
            frames.append(df)
        return pd.concat(frames, ignore_index=True)

    # ------------------------------------------------------------------ #
    #  Track plots / contacts / logs                                       #
    # ------------------------------------------------------------------ #

    def track_plots(
        self,
        *,
        method: _MethodArg | None = None,
        scale: _ScaleArg | None = None,
        spikein_method: _SpikeInArg | None = None,
        merged: bool | None = None,
        only_existing: bool = True,
    ) -> list[Path]:
        """Return PlotNado track plot files.

        Track plots live under ``track_plots/{method}/{scale}/`` or
        ``track_plots/merged/{method}/{scale}/`` and may use any image format
        (svg, png, pdf) depending on the pipeline configuration.

        ``method``, ``scale``, and ``spikein_method`` accept enums or strings.
        """
        tp_dir = self._output_dir / "track_plots"
        if not tp_dir.exists():
            return []

        paths = list(tp_dir.rglob("*"))
        paths = [p for p in paths if p.is_file() and p.suffix in {".svg", ".png", ".pdf"}]

        method_s = _ev(method)
        scale_s = _ev(scale)
        spikein_s = _ev(spikein_method)

        if method_s is not None:
            paths = [p for p in paths if f"/{method_s}/" in p.as_posix()]
        if spikein_s is not None:
            paths = [p for p in paths if f"/spikein/{spikein_s}/" in p.as_posix()]
        elif scale_s is not None:
            paths = [p for p in paths if f"/{scale_s}/" in p.as_posix()]
        if merged is True:
            paths = [p for p in paths if "/merged/" in p.as_posix()]
        elif merged is False:
            paths = [p for p in paths if "/merged/" not in p.as_posix()]

        if only_existing:
            paths = [p for p in paths if p.exists()]
        return sorted(paths)

    def contacts(self, only_existing: bool = True) -> list[Path]:
        """Return MCC contact matrix paths (``.mcool`` files)."""
        contacts_dir = self._output_dir / "mcc" / "contacts"
        if not contacts_dir.exists():
            return []
        paths = sorted(contacts_dir.rglob("*.mcool"))
        if only_existing:
            paths = [p for p in paths if p.exists()]
        return paths

    def logs(
        self,
        rule: str | None = None,
        sample: str | None = None,
        only_existing: bool = True,
    ) -> list[Path]:
        """Return Snakemake log file paths.

        Parameters
        ----------
        rule:   Rule name directory, e.g. ``"bowtie2"``, ``"macs2_callpeak"``.
        sample: Sample name stem to filter within the rule directory.
        """
        log_dir = self._output_dir / "logs"
        if not log_dir.exists():
            return []

        if rule is not None:
            search_dir = log_dir / rule
            if not search_dir.exists():
                return []
            paths = sorted(search_dir.glob("*.log"))
        else:
            paths = sorted(log_dir.rglob("*.log"))

        if sample is not None:
            paths = [p for p in paths if p.stem == sample or p.stem.startswith(f"{sample}.")]

        if only_existing:
            paths = [p for p in paths if p.exists()]
        return paths

    # ------------------------------------------------------------------ #
    #  Load QC / performance data                                          #
    # ------------------------------------------------------------------ #

    def load_alignment_stats(self) -> pd.DataFrame:
        """Load per-sample alignment step counts as a ``pd.DataFrame``.

        Reads ``qc/alignment_stats.tsv`` which is produced by the
        ``alignment_stats`` aggregation rule.  Columns include ``sample``
        plus one column per processing step (reads lost / retained).
        """
        p = self._output_dir / "qc" / "alignment_stats.tsv"
        if not p.exists():
            raise FileNotFoundError(f"alignment_stats.tsv not found at {p}")
        return pd.read_csv(p, sep="\t")

    def load_frip(self) -> pd.DataFrame:
        """Load FRiP (Fraction of Reads in Peaks) scores as a ``pd.DataFrame``.

        Reads all ``qc/frip_enrichment/{method}/{sample}_frip.txt`` files and
        concatenates them with ``sample`` and ``peak_method`` columns added.
        The files are deeptools ``plotEnrichment --outRawCounts`` output.
        """
        frip_dir = self._output_dir / "qc" / "frip_enrichment"
        if not frip_dir.exists():
            raise FileNotFoundError(f"FRiP directory not found: {frip_dir}")

        frames = []
        for p in sorted(frip_dir.rglob("*_frip.txt")):
            try:
                df = pd.read_csv(p, sep="\t", comment="#")
                # Strip leading '#' from header if present
                df.columns = [c.lstrip("#").strip() for c in df.columns]
                # Derive sample and peak_method from path:
                # frip_enrichment/{peak_method}/{sample}_frip.txt
                peak_method = p.parent.name
                sample_name = p.stem.removesuffix("_frip")
                df.insert(0, "sample", sample_name)
                df.insert(1, "peak_method", peak_method)
                frames.append(df)
            except Exception as e:
                import warnings
                warnings.warn(f"Failed to parse FRiP file {p}: {e}", stacklevel=2)
                continue

        if not frames:
            raise FileNotFoundError("No FRiP files found or parseable.")
        return pd.concat(frames, ignore_index=True)

    def load_library_complexity(self) -> pd.DataFrame:
        """Load Picard MarkDuplicates library complexity metrics as a ``pd.DataFrame``.

        Reads all ``qc/library_complexity/{sample}.metrics`` files.
        Returns the METRICS section rows with a ``sample`` column prepended.
        """
        lc_dir = self._output_dir / "qc" / "library_complexity"
        if not lc_dir.exists():
            raise FileNotFoundError(f"Library complexity directory not found: {lc_dir}")

        _SAMTOOLS_HEADER = {"FILENAME", "COMMAND", "FRAMES", "READ_EXAMINED", "DUPLICATE_TOTAL"}

        frames = []
        for p in sorted(lc_dir.glob("*.metrics")):
            sample_name = p.stem
            try:
                import io
                lines = p.read_text().splitlines()
                non_empty = [ln for ln in lines if ln.strip()]

                # samtools markdup TSV: first line is tab-separated header
                if non_empty and _SAMTOOLS_HEADER & set(non_empty[0].split("\t")):
                    df = pd.read_csv(io.StringIO("\n".join(non_empty)), sep="\t")
                    df = df.dropna(how="all")
                else:
                    # Picard metrics files: header block then TSV starting with "LIBRARY"
                    header_idx = next(
                        (i for i, line in enumerate(lines)
                         if line.startswith("LIBRARY") or line.startswith("ESTIMATED_LIBRARY")),
                        None,
                    )
                    if header_idx is not None:
                        df = pd.read_csv(io.StringIO("\n".join(lines[header_idx:])), sep="\t")
                        df = df.dropna(how="all")
                    else:
                        # Fallback: skip comment/blank lines
                        data_lines = [ln for ln in lines if ln and not ln.startswith("#")]
                        if len(data_lines) < 2:
                            continue
                        df = pd.read_csv(io.StringIO("\n".join(data_lines)), sep="\t")

                df.insert(0, "sample", sample_name)
                frames.append(df)
            except Exception:
                continue

        if not frames:
            raise FileNotFoundError("No library complexity metrics files found or parseable.")
        return pd.concat(frames, ignore_index=True)

    def load_benchmarks(self, rule: str | None = None) -> pd.DataFrame:
        """Load Snakemake benchmark timing / memory data as a ``pd.DataFrame``.

        Reads all ``.benchmark/{rule}/{sample}.tsv`` files.  Adds ``rule``
        and ``sample`` columns derived from the path.

        Parameters
        ----------
        rule: Limit to a single rule directory, e.g. ``"bowtie2"``.
        """
        bench_dir = self._output_dir / ".benchmark"
        if not bench_dir.exists():
            raise FileNotFoundError(f"Benchmark directory not found: {bench_dir}")

        if rule is not None:
            search_dir = bench_dir / rule
            tsv_files = sorted(search_dir.rglob("*.tsv")) if search_dir.exists() else []
        else:
            tsv_files = sorted(bench_dir.rglob("*.tsv"))

        frames = []
        for p in tsv_files:
            try:
                df = pd.read_csv(p, sep="\t")
                df.insert(0, "rule", p.parent.name)
                df.insert(1, "sample", p.stem)
                frames.append(df)
            except Exception:
                continue

        if not frames:
            raise FileNotFoundError("No benchmark files found or parseable.")
        return pd.concat(frames, ignore_index=True)

    def load_peaks(
        self,
        *,
        sample: str | list[str] | None = None,
        condition: str | list[str] | None = None,
        antibody: str | list[str] | None = None,
        group: str | list[str] | None = None,
        method: str | None = None,
        merged: bool | None = None,
        extra_cols: list[str] | None = None,
    ) -> pd.DataFrame:
        """Load peak BED files as a concatenated ``pd.DataFrame``.

        Applies the same filters as :meth:`peaks`.  Adds ``sample``,
        ``peak_method``, and ``merged`` columns.  BED columns are named
        ``chrom``, ``start``, ``end``, ``name``, ``score``, ``strand``
        (extra columns numbered ``col_7``, ``col_8``, …).

        Parameters
        ----------
        extra_cols: Override names for columns beyond the standard 6.
        """
        paths = self.peaks(
            sample=sample, condition=condition, antibody=antibody,
            group=group, method=method, merged=merged,
        )
        if not paths:
            return pd.DataFrame()

        peaks_dir = self._output_dir / "peaks"
        std_cols = ["chrom", "start", "end", "name", "score", "strand"]
        frames = []

        for p in paths:
            try:
                df = pd.read_csv(p, sep="\t", header=None, comment="#")
                n = len(df.columns)
                if n <= len(std_cols):
                    df.columns = std_cols[:n]
                else:
                    extra = extra_cols or [f"col_{i+7}" for i in range(n - len(std_cols))]
                    df.columns = std_cols + extra[: n - len(std_cols)]
                # Derive metadata from path
                rec = _parse_peak_path(p, peaks_dir) or {}
                df["sample"] = rec.get("sample", p.stem)
                df["peak_method"] = rec.get("method", "unknown")
                df["merged"] = rec.get("merged", False)
                frames.append(df)
            except Exception:
                continue

        return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()

    # ------------------------------------------------------------------ #
    #  Convenience helpers                                                 #
    # ------------------------------------------------------------------ #

    def protocol(self) -> Path | None:
        """Path to the auto-generated ``protocol.txt`` file."""
        p = self._output_dir / "protocol.txt"
        return p if p.exists() else None

    def bai_for(self, bam: str | Path) -> Path:
        """Return the BAI index path for a BAM file.

        Checks both ``{bam}.bai`` and ``{bam_stem}.bai`` conventions.
        Raises ``FileNotFoundError`` if neither exists.
        """
        bam = Path(bam)
        candidates = [Path(str(bam) + ".bai"), bam.with_suffix(".bai")]
        for c in candidates:
            if c.exists():
                return c
        raise FileNotFoundError(f"No BAI index found for {bam}")

    def sample_files(self, sample: str) -> dict:
        """Return a dict of all existing files for a single sample.

        Keys: ``bam``, ``bai``, ``bigwigs``, ``peaks``, ``qc``,
        ``methylation``, ``vcf``.
        """
        bam_list = self.bams(sample=sample)
        bam = bam_list[0] if bam_list else None
        try:
            bai = self.bai_for(bam) if bam else None
        except FileNotFoundError:
            bai = None

        return {
            "bam": bam,
            "bai": bai,
            "bigwigs": self.bigwigs(sample=sample),
            "peaks": self.peaks(sample=sample),
            "qc": self.qc(sample=sample),
            "methylation": self.methylation(sample=sample),
            "vcf": self.vcf(sample=sample),
        }

    def condition_pairs(self) -> list[tuple[str, str]]:
        """Return all ordered condition pairs for pairwise comparisons.

        Useful when constructing DESeq2 contrasts or subtraction bigWig paths.

        Returns
        -------
        list of ``(condition_a, condition_b)`` tuples (all permutations).
        """
        from itertools import permutations
        return list(permutations(self.conditions, 2))

    # ------------------------------------------------------------------ #
    #  Filtering view                                                      #
    # ------------------------------------------------------------------ #

    def filter(
        self,
        *,
        condition: str | list[str] | None = None,
        antibody: str | list[str] | None = None,
        sample: str | list[str] | None = None,
        group: str | list[str] | None = None,
    ) -> "_FilteredProject":
        """Return a filtered view of the project.

        All file-access methods on the returned object are pre-filtered
        to the matching samples.

        Example
        -------
        ::

            treated = proj.filter(condition="treated")
            treated.bigwigs(method="deeptools")
            treated.peaks(method="macs2")
        """
        targets = self._resolve_samples(sample, condition, antibody, group)
        return _FilteredProject(self, targets)

    # ------------------------------------------------------------------ #
    #  Diagnostics                                                         #
    # ------------------------------------------------------------------ #

    # ------------------------------------------------------------------ #
    #  Escape hatches / utilities                                         #
    # ------------------------------------------------------------------ #

    def reload(self) -> None:
        """Clear all cached properties so they are recomputed on next access."""
        for name in (
            "_bw_index", "_peak_index", "_sample_meta",
            "samples", "conditions", "antibodies", "groups", "assay",
            "pileup_methods", "scales", "spikein_methods", "peak_methods",
            "consensus_groups", "log_rules", "benchmark_rules",
        ):
            self.__dict__.pop(name, None)

    def select_files(self, pattern: str, recursive: bool = True) -> list[Path]:
        """Glob for files matching *pattern* under the output directory.

        Parameters
        ----------
        pattern:   Glob pattern, e.g. ``"logs/bowtie2/*.log"`` or ``"**/*.bam"``.
        recursive: When ``True`` (default) and the pattern contains no ``**``,
                   use ``rglob``; otherwise ``glob`` is used as-is.
        """
        if "**" in pattern or not recursive:
            return sorted(self._output_dir.glob(pattern))
        return sorted(self._output_dir.rglob(pattern))

    def enrich(self, paths: list[Path]) -> pd.DataFrame:
        """Annotate file paths with sample metadata from the design.

        Extracts sample names from file stems (stripping ``_plus`` / ``_minus``
        strand suffixes) and left-joins ``condition``, ``antibody``, ``group``
        from the design.

        Parameters
        ----------
        paths: Iterable of ``pathlib.Path`` objects (e.g. from :meth:`bigwigs`).
        """
        records = []
        for p in paths:
            stem = p.stem
            for suffix in ("_plus", "_minus"):
                if stem.endswith(suffix):
                    stem = stem[: -len(suffix)]
                    break
            records.append({"path": p, "sample": stem})

        df = pd.DataFrame(records)
        if df.empty:
            return df

        meta = self._sample_meta.reset_index()
        return df.merge(meta, on="sample", how="left")

    # ------------------------------------------------------------------ #
    #  DE / CRISPR loaders                                                #
    # ------------------------------------------------------------------ #

    def load_deseq2(self, contrast: str | None = None) -> pd.DataFrame:
        """Load DESeq2 differential expression results as a ``pd.DataFrame``.

        Searches ``{output_dir}/deseq2_results/DEseq2_*.csv`` and the parent
        directory (where Quarto renders by default).  Normalised-counts CSVs
        are excluded.

        Parameters
        ----------
        contrast: Load a single named contrast.  When ``None``, all contrasts
                  are concatenated with a ``contrast`` column added.
        """
        search_dirs = [
            self._output_dir / "deseq2_results",
            self._output_dir.parent,
        ]
        if contrast is not None:
            for d in search_dirs:
                p = d / f"DEseq2_{contrast}.csv"
                if p.exists():
                    df = pd.read_csv(p)
                    df.insert(0, "contrast", contrast)
                    return df
            raise FileNotFoundError(
                f"DESeq2 results for contrast '{contrast}' not found in {search_dirs}."
            )

        frames = []
        seen: set[str] = set()
        for d in search_dirs:
            if not d.exists():
                continue
            for p in sorted(d.glob("DEseq2_*.csv")):
                name = p.stem.removeprefix("DEseq2_")
                if "_normalised_counts" in name or p.stem in seen:
                    continue
                seen.add(p.stem)
                df = pd.read_csv(p)
                df.insert(0, "contrast", name)
                frames.append(df)

        if not frames:
            raise FileNotFoundError("No DESeq2 result CSVs found.")
        return pd.concat(frames, ignore_index=True)

    def load_mageck(
        self,
        type: str = "gene",
        analysis: str = "test",
    ) -> pd.DataFrame:
        """Load a MAGeCK CRISPR screen summary as a ``pd.DataFrame``.

        Reads ``readcounts/mageck/mageck_{analysis}.{type}_summary.txt``.

        Parameters
        ----------
        type:     ``"gene"`` (default) or ``"sgrna"``.
        analysis: ``"test"`` (default) or ``"mle"``.
        """
        p = (
            self._output_dir
            / "readcounts"
            / "mageck"
            / f"mageck_{analysis}.{type}_summary.txt"
        )
        if not p.exists():
            raise FileNotFoundError(f"MAGeCK summary not found: {p}")
        return pd.read_csv(p, sep="\t")

    # ------------------------------------------------------------------ #
    #  Diagnostics                                                         #
    # ------------------------------------------------------------------ #

    def summary(self, detail: bool = False) -> pd.DataFrame:
        """Return a DataFrame summarising which file types are present.

        Parameters
        ----------
        detail: When ``True``, expand BigWig and peak rows to show
                per-(method, scale) and per-method counts.
        """
        if detail:
            rows: list[dict] = []
            bam_count = len(self.bams())
            rows.append({
                "File type": "BAMs", "Detail": "all samples",
                "Count": bam_count, "Present": "✓" if bam_count else "–",
            })

            if not self._bw_index.empty:
                for (method, scale), grp in self._bw_index.groupby(["method", "scale"]):
                    rows.append({
                        "File type": "BigWig", "Detail": f"{method} / {scale}",
                        "Count": len(grp), "Present": "✓",
                    })
            else:
                rows.append({"File type": "BigWig", "Detail": "–", "Count": 0, "Present": "–"})

            if not self._peak_index.empty:
                for method, grp in self._peak_index.groupby("method"):
                    rows.append({
                        "File type": "Peaks", "Detail": method,
                        "Count": len(grp), "Present": "✓",
                    })
            else:
                rows.append({"File type": "Peaks", "Detail": "–", "Count": 0, "Present": "–"})

            for name, present in [
                ("Counts",         self.counts() is not None),
                ("QC reports",     bool(self.qc())),
                ("MultiQC report", self.report() is not None),
                ("Methylation",    bool(self.methylation())),
                ("VCF",            bool(self.vcf())),
                ("Heatmaps",       bool(self.heatmaps())),
                ("Norm. factors",  bool(self.normalisation_factors())),
            ]:
                rows.append({
                    "File type": name, "Detail": "–",
                    "Count": None, "Present": "✓" if present else "–",
                })
            return pd.DataFrame(rows)

        rows_simple = []
        checks = [
            ("BAMs",            bool(self.bams())),
            ("BigWigs",         not self._bw_index.empty),
            ("Peaks",           not self._peak_index.empty),
            ("Counts",          self.counts() is not None),
            ("QC reports",      bool(self.qc())),
            ("MultiQC report",  self.report() is not None),
            ("Methylation",     bool(self.methylation())),
            ("VCF",             bool(self.vcf())),
            ("Heatmaps",        bool(self.heatmaps())),
            ("Norm. factors",   bool(self.normalisation_factors())),
        ]
        for name, present in checks:
            rows_simple.append({"File type": name, "Present": "✓" if present else "–"})
        return pd.DataFrame(rows_simple)

    def __repr__(self) -> str:
        return (
            f"SeqNadoProject("
            f"assay={self.assay!r}, "
            f"output_dir={str(self._output_dir)!r}, "
            f"samples={len(self.samples)}, "
            f"conditions={self.conditions}"
            f")"
        )


# ---------------------------------------------------------------------------
# Filtered view
# ---------------------------------------------------------------------------

class _FilteredProject:
    """Project view pre-filtered to a subset of samples.

    Created by :meth:`SeqNadoProject.filter`.
    """

    def __init__(self, project: SeqNadoProject, targets: set[str] | None):
        self._project = project
        self._targets = targets

    def _sample_kwarg(self) -> list[str] | None:
        return list(self._targets) if self._targets is not None else None

    @property
    def samples(self) -> list[str]:
        if self._targets is None:
            return self._project.samples
        return sorted(self._targets & set(self._project.samples))

    @property
    def conditions(self) -> list[str]:
        meta = self._project._sample_meta
        if self._targets is None or meta.empty:
            return self._project.conditions
        sub = meta[meta.index.isin(self._targets)]
        return sorted(sub["condition"].dropna().unique().tolist())

    @property
    def antibodies(self) -> list[str]:
        meta = self._project._sample_meta
        if self._targets is None or meta.empty:
            return self._project.antibodies
        sub = meta[meta.index.isin(self._targets)]
        return sorted(sub["antibody"].dropna().unique().tolist())

    @property
    def pileup_methods(self) -> list[str]:
        df = self._project._bw_index
        if df.empty or self._targets is None:
            return self._project.pileup_methods
        return sorted(df[df["sample"].isin(self._targets)]["method"].dropna().unique().tolist())

    @property
    def scales(self) -> list[str]:
        df = self._project._bw_index
        if df.empty or self._targets is None:
            return self._project.scales
        return sorted(df[df["sample"].isin(self._targets)]["scale"].dropna().unique().tolist())

    @property
    def peak_methods(self) -> list[str]:
        df = self._project._peak_index
        if df.empty or self._targets is None:
            return self._project.peak_methods
        return sorted(df[df["sample"].isin(self._targets)]["method"].dropna().unique().tolist())

    @property
    def consensus_groups(self) -> list[str]:
        return self._project.consensus_groups

    @property
    def design(self) -> pd.DataFrame | None:
        df = self._project.design
        if df is None or self._targets is None:
            return df
        if "sample_id" in df.columns:
            return df[df["sample_id"].isin(self._targets)].copy()
        return df

    def bigwigs(self, **kwargs) -> list[Path]:
        kwargs.setdefault("sample", self._sample_kwarg())
        return self._project.bigwigs(**kwargs)

    def peaks(self, **kwargs) -> list[Path]:
        kwargs.setdefault("sample", self._sample_kwarg())
        return self._project.peaks(**kwargs)

    def bams(self, **kwargs) -> list[Path]:
        kwargs.setdefault("sample", self._sample_kwarg())
        return self._project.bams(**kwargs)

    def methylation(self, **kwargs) -> list[Path]:
        kwargs.setdefault("sample", self._sample_kwarg())
        return self._project.methylation(**kwargs)

    def vcf(self, **kwargs) -> list[Path]:
        kwargs.setdefault("sample", self._sample_kwarg())
        return self._project.vcf(**kwargs)

    def qc(self, **kwargs) -> list[Path]:
        return self._project.qc(**kwargs)

    def track_plots(self, **kwargs) -> list[Path]:
        return self._project.track_plots(**kwargs)

    def contacts(self, **kwargs) -> list[Path]:
        return self._project.contacts(**kwargs)

    def logs(self, **kwargs) -> list[Path]:
        return self._project.logs(**kwargs)

    def load_peaks(self, **kwargs) -> "pd.DataFrame":
        kwargs.setdefault("sample", self._sample_kwarg())
        return self._project.load_peaks(**kwargs)

    def sample_files(self, sample: str) -> dict:
        return self._project.sample_files(sample)

    def condition_pairs(self) -> list[tuple[str, str]]:
        return self._project.condition_pairs()

    def filter(
        self,
        *,
        condition: str | list[str] | None = None,
        antibody: str | list[str] | None = None,
        sample: str | list[str] | None = None,
        group: str | list[str] | None = None,
    ) -> "_FilteredProject":
        extra = self._project._resolve_samples(sample, condition, antibody, group)
        if self._targets is None:
            return _FilteredProject(self._project, extra)
        merged = self._targets if extra is None else self._targets & extra
        return _FilteredProject(self._project, merged)

    def __repr__(self) -> str:
        n = len(self._targets) if self._targets is not None else len(self._project.samples)
        return f"_FilteredProject(samples={n}, from={self._project!r})"


# ---------------------------------------------------------------------------
# Multiomics wrapper
# ---------------------------------------------------------------------------

def _is_multiomics_root(path: Path) -> bool:
    """Return True if *path* looks like a multiomics seqnado_output root.

    Detected when it contains at least one subdirectory whose name matches
    a known assay clean-name (atac, chip, rna, …) AND that subdirectory
    itself contains an ``aligned/`` directory.
    """
    from seqnado.core import Assay
    known = set(Assay.all_assay_clean_names())
    for child in path.iterdir():
        if child.is_dir() and child.name in known:
            if (child / "aligned").exists() or (child / "bigwigs").exists():
                return True
    return False


class SeqNadoMultiProject:
    """Analysis interface for a multiomics SeqNado output directory.

    Multiomics runs place each assay's outputs under
    ``seqnado_output/{assay}/``, e.g.:

    .. code-block:: text

        seqnado_output/
        ├── chip/
        ├── atac/
        └── rna/

    ``SeqNadoMultiProject`` wraps one :class:`SeqNadoProject` per detected
    assay and exposes them via attribute / item access.

    Parameters
    ----------
    output_dir:
        Top-level output directory (e.g. ``seqnado_output/``).

    Examples
    --------
    ::

        mp = SeqNadoMultiProject("seqnado_output/")
        mp.assays              # ['atac', 'chip', 'rna']
        mp.chip                # SeqNadoProject for chip sub-directory
        mp["atac"]             # SeqNadoProject for atac sub-directory

        mp.chip.bigwigs(scale="csaw")
        mp.rna.load_counts()

        # Iterate all assays
        for name, proj in mp.items():
            print(name, proj.samples)
    """

    def __init__(self, output_dir: str | Path):
        self._root = Path(output_dir)
        from seqnado.core import Assay
        known = set(Assay.all_assay_clean_names())
        self._projects: dict[str, SeqNadoProject] = {}
        for child in sorted(self._root.iterdir()):
            if child.is_dir() and child.name in known:
                if (child / "aligned").exists() or (child / "bigwigs").exists():
                    self._projects[child.name] = SeqNadoProject(child)

        if not self._projects:
            raise ValueError(
                f"{self._root} does not appear to be a multiomics output root "
                f"(no assay subdirectories with aligned/ or bigwigs/ found)."
            )

    @property
    def assays(self) -> list[str]:
        """Assay names present in the output directory."""
        return list(self._projects)

    @property
    def samples(self) -> list[str]:
        """All sample names across all assays (union, sorted)."""
        s: set[str] = set()
        for proj in self._projects.values():
            s.update(proj.samples)
        return sorted(s)

    @property
    def conditions(self) -> list[str]:
        """All condition labels across all assays (union, sorted)."""
        s: set[str] = set()
        for proj in self._projects.values():
            s.update(proj.conditions)
        return sorted(s)

    @property
    def antibodies(self) -> list[str]:
        """All antibody / IP targets across all assays (union, sorted)."""
        s: set[str] = set()
        for proj in self._projects.values():
            s.update(proj.antibodies)
        return sorted(s)

    @property
    def groups(self) -> list[str]:
        """All scaling groups across all assays (union, sorted)."""
        s: set[str] = set()
        for proj in self._projects.values():
            s.update(proj.groups)
        return sorted(s)

    def __getitem__(self, assay: str) -> SeqNadoProject:
        if assay not in self._projects:
            raise KeyError(f"Assay '{assay}' not found. Available: {self.assays}")
        return self._projects[assay]

    def __getattr__(self, name: str) -> SeqNadoProject:
        if name.startswith("_"):
            raise AttributeError(name)
        if name in self._projects:
            return self._projects[name]
        raise AttributeError(
            f"No assay '{name}' found. Available: {self.assays}"
        )

    def items(self):
        """Iterate ``(assay_name, SeqNadoProject)`` pairs."""
        return self._projects.items()

    def values(self):
        return self._projects.values()

    # ------------------------------------------------------------------ #
    #  Cross-assay file access                                             #
    # ------------------------------------------------------------------ #

    def _cross(
        self,
        method_name: str,
        assays: list[str] | None,
        kwargs: dict,
    ) -> pd.DataFrame:
        """Call *method_name* on each per-assay project and concatenate results.

        Returns a DataFrame with an ``assay`` column prepended.  File-access
        methods that return ``list[Path]`` are expanded into a ``path`` column;
        methods that already return a DataFrame have the column inserted.
        """
        target_assays = assays if assays is not None else self.assays
        frames = []
        for name in target_assays:
            proj = self._projects.get(name)
            if proj is None:
                continue
            result = getattr(proj, method_name)(**kwargs)
            if isinstance(result, list):
                if not result:
                    continue
                df = pd.DataFrame({"path": result})
                df.insert(0, "assay", name)
            else:
                if result.empty:
                    continue
                df = result.copy()
                df.insert(0, "assay", name)
            frames.append(df)
        return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()

    def bigwigs(
        self,
        *,
        assays: list[str] | None = None,
        sample: str | list[str] | None = None,
        condition: str | list[str] | None = None,
        antibody: str | list[str] | None = None,
        group: str | list[str] | None = None,
        method: "_MethodArg | None" = None,
        scale: "_ScaleArg | None" = None,
        spikein_method: "_SpikeInArg | None" = None,
        merged: bool | None = None,
        strand: str | None = None,
        only_existing: bool = True,
    ) -> pd.DataFrame:
        """Return bigWig paths across all assays as a DataFrame.

        Accepts the same filter kwargs as :meth:`SeqNadoProject.bigwigs` plus
        an optional ``assays`` list to restrict which assays are queried.

        Returns a ``pd.DataFrame`` with an ``assay`` column followed by a
        ``path`` column.

        Examples
        --------
        ::

            mp.bigwigs(condition="DMSO")
            mp.bigwigs(condition="DMSO", assays=["chip", "cat"])
            mp.bigwigs(scale="unscaled", merged=False)
        """
        kw = dict(
            sample=sample, condition=condition, antibody=antibody, group=group,
            method=method, scale=scale, spikein_method=spikein_method,
            merged=merged, strand=strand, only_existing=only_existing,
        )
        return self._cross("bigwigs", assays, kw)

    def peaks(
        self,
        *,
        assays: list[str] | None = None,
        sample: str | list[str] | None = None,
        condition: str | list[str] | None = None,
        antibody: str | list[str] | None = None,
        group: str | list[str] | None = None,
        method: "_PeakMethodArg | None" = None,
        merged: bool | None = None,
        only_existing: bool = True,
    ) -> pd.DataFrame:
        """Return peak paths across all assays as a DataFrame.

        Same filters as :meth:`SeqNadoProject.peaks` plus optional ``assays``.
        Returns ``pd.DataFrame`` with ``assay`` and ``path`` columns.
        """
        kw = dict(
            sample=sample, condition=condition, antibody=antibody, group=group,
            method=method, merged=merged, only_existing=only_existing,
        )
        return self._cross("peaks", assays, kw)

    def bams(
        self,
        *,
        assays: list[str] | None = None,
        sample: str | list[str] | None = None,
        condition: str | list[str] | None = None,
        antibody: str | list[str] | None = None,
        group: str | list[str] | None = None,
        merged: bool = False,
        only_existing: bool = True,
    ) -> pd.DataFrame:
        """Return BAM paths across all assays as a DataFrame.

        Same filters as :meth:`SeqNadoProject.bams` plus optional ``assays``.
        Returns ``pd.DataFrame`` with ``assay`` and ``path`` columns.
        """
        kw = dict(
            sample=sample, condition=condition, antibody=antibody, group=group,
            merged=merged, only_existing=only_existing,
        )
        return self._cross("bams", assays, kw)

    def load_peaks(
        self,
        *,
        assays: list[str] | None = None,
        **kwargs,
    ) -> pd.DataFrame:
        """Load and concatenate peak BED files across assays.

        Returns a DataFrame with an ``assay`` column plus the standard peak
        columns (``chrom``, ``start``, ``end``, ``sample``, ``peak_method``, …).
        """
        return self._cross("load_peaks", assays, kwargs)

    # ------------------------------------------------------------------ #
    #  Cross-assay filter                                                  #
    # ------------------------------------------------------------------ #

    def filter(
        self,
        *,
        condition: str | list[str] | None = None,
        antibody: str | list[str] | None = None,
        sample: str | list[str] | None = None,
        group: str | list[str] | None = None,
        assays: list[str] | None = None,
    ) -> "_MultiFilteredProject":
        """Return a filtered cross-assay view.

        The returned :class:`_MultiFilteredProject` pre-applies sample /
        condition / antibody filters to every file-access method **and**
        optionally restricts which assays are included.

        Examples
        --------
        ::

            mp.filter(condition="DMSO").bigwigs()
            mp.filter(condition="DMSO", assays=["chip", "cat"]).peaks()
        """
        target_assays = assays if assays is not None else self.assays
        filtered: dict[str, "_FilteredProject"] = {}
        for name in target_assays:
            proj = self._projects.get(name)
            if proj is not None:
                filtered[name] = proj.filter(
                    condition=condition, antibody=antibody,
                    sample=sample, group=group,
                )
        return _MultiFilteredProject(self, filtered)

    # ------------------------------------------------------------------ #
    #  Diagnostics                                                         #
    # ------------------------------------------------------------------ #

    def summary(self) -> pd.DataFrame:
        """Combined summary table across all assays."""
        frames = []
        for name, proj in self._projects.items():
            df = proj.summary()
            df.insert(0, "assay", name)
            frames.append(df)
        return pd.concat(frames, ignore_index=True)

    def what_exists(self) -> pd.DataFrame:
        """Per-sample availability table across all assays."""
        frames = []
        for name, proj in self._projects.items():
            df = proj.what_exists()
            df.insert(0, "assay", name)
            frames.append(df)
        return pd.concat(frames, ignore_index=True)

    def __repr__(self) -> str:
        return f"SeqNadoMultiProject(root={str(self._root)!r}, assays={self.assays})"


# ---------------------------------------------------------------------------
# Cross-assay filtered view
# ---------------------------------------------------------------------------

class _MultiFilteredProject:
    """Cross-assay view pre-filtered to a subset of samples.

    Created by :meth:`SeqNadoMultiProject.filter`.
    """

    def __init__(
        self,
        multi: SeqNadoMultiProject,
        filtered: "dict[str, _FilteredProject]",
    ):
        self._multi = multi
        self._filtered = filtered

    @property
    def assays(self) -> list[str]:
        return list(self._filtered)

    def _cross(self, method_name: str, kwargs: dict) -> pd.DataFrame:
        frames = []
        for name, fproj in self._filtered.items():
            result = getattr(fproj, method_name)(**kwargs)
            if isinstance(result, list):
                if not result:
                    continue
                df = pd.DataFrame({"path": result})
                df.insert(0, "assay", name)
            else:
                if result.empty:
                    continue
                df = result.copy()
                df.insert(0, "assay", name)
            frames.append(df)
        return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()

    def bigwigs(self, **kwargs) -> pd.DataFrame:
        return self._cross("bigwigs", kwargs)

    def peaks(self, **kwargs) -> pd.DataFrame:
        return self._cross("peaks", kwargs)

    def bams(self, **kwargs) -> pd.DataFrame:
        return self._cross("bams", kwargs)

    def load_peaks(self, **kwargs) -> pd.DataFrame:
        return self._cross("load_peaks", kwargs)

    def filter(
        self,
        *,
        condition: str | list[str] | None = None,
        antibody: str | list[str] | None = None,
        sample: str | list[str] | None = None,
        group: str | list[str] | None = None,
        assays: list[str] | None = None,
    ) -> "_MultiFilteredProject":
        target = assays if assays is not None else self.assays
        new_filtered = {}
        for name in target:
            fproj = self._filtered.get(name)
            if fproj is not None:
                new_filtered[name] = fproj.filter(
                    condition=condition, antibody=antibody,
                    sample=sample, group=group,
                )
        return _MultiFilteredProject(self._multi, new_filtered)

    def __repr__(self) -> str:
        counts = {n: len(f.samples) for n, f in self._filtered.items()}
        return f"_MultiFilteredProject(assays={counts})"


def open_project(output_dir: str | Path, **kwargs) -> "SeqNadoProject | SeqNadoMultiProject":
    """Open a SeqNado output directory, auto-detecting single vs multiomics.

    Returns a :class:`SeqNadoMultiProject` when the directory is a
    multiomics root (contains assay subdirectories), otherwise a
    :class:`SeqNadoProject`.

    Parameters
    ----------
    output_dir: Path to the output directory.
    **kwargs:   Forwarded to :class:`SeqNadoProject` (``config``, ``design``).
    """
    path = Path(output_dir)
    if _is_multiomics_root(path):
        return SeqNadoMultiProject(path)
    return SeqNadoProject(path, **kwargs)
