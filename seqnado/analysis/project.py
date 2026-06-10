"""Post-pipeline analysis API for SeqNado output directories.

Usage:

    from seqnado.analysis import SeqNadoProject

    proj = SeqNadoProject("seqnado_output")
    # or with explicit config
    proj = SeqNadoProject("seqnado_output", config="chip_config.yml")

    proj.samples          # list[str]
    proj.conditions       # list[str]
    proj.antibodies       # list[str]  (ChIP / CUT&Tag)
    proj.design           # pd.DataFrame

    proj.bigwigs()                                      # all bigWigs
    proj.bigwigs(condition="treated")                   # by condition
    proj.bigwigs(antibody="H3K27ac", scale="unscaled")  # by antibody + scale
    proj.bigwigs(method="deeptools", merged=True)       # merged tracks

    proj.peaks(method="macs2")
    proj.bams()
    proj.counts()                    # Path to counts TSV
    proj.load_counts()               # pd.DataFrame

    proj.filter(condition="treated").bigwigs()
    proj.filter(antibody="H3K27ac").peaks(merged=True)

    proj.bigwig_dataframe()          # DataFrame with path + metadata columns
    proj.peak_dataframe()
"""

from __future__ import annotations

from functools import cached_property
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd

if TYPE_CHECKING:
    from seqnado.config import SeqnadoConfig


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
            return config

        path: Path | None = None
        if config is not None:
            path = Path(config)
        else:
            # auto-discover in project dir (parent of output_dir) or cwd
            for search_dir in [self._output_dir.parent, Path(".")]:
                candidates = sorted(search_dir.glob("*_config.yml"))
                if candidates:
                    path = candidates[0]
                    break

        if path is None or not path.exists():
            return None

        try:
            data = _load_yaml(path)
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
            path = Path(self._config.metadata)
        else:
            for candidate in [
                self._output_dir.parent / "design.tsv",
                self._output_dir.parent / "design.csv",
                self._output_dir / "design.tsv",
                self._output_dir / "design.csv",
            ]:
                if candidate.exists():
                    path = candidate
                    break

        if path is None or not path.exists():
            return None

        sep = "," if path.suffix == ".csv" else "\t"
        return pd.read_csv(path, sep=sep)

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
        if self._design_df is None:
            return pd.DataFrame(columns=["condition", "antibody", "group"])

        df = self._design_df.copy()

        # Normalise column names
        if "ip" in df.columns and "antibody" not in df.columns:
            df = df.rename(columns={"ip": "antibody"})

        rows = []
        for _, row in df.iterrows():
            sid = str(row["sample_id"])
            cond = row.get("condition") if "condition" in df.columns else None
            ab = row.get("antibody") if "antibody" in df.columns else None
            grp = row.get("scaling_group") if "scaling_group" in df.columns else None

            # Non-IP sample (ATAC, RNA, etc.)
            rows.append({"sample": sid, "condition": cond, "antibody": ab, "group": grp})

            # IP-based sample: full name is sample_id + "_" + antibody
            if pd.notna(ab) and ab:
                rows.append({"sample": f"{sid}_{ab}", "condition": cond, "antibody": ab, "group": grp})

        meta = pd.DataFrame(rows).drop_duplicates("sample").set_index("sample")
        return meta

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
        method: str | None = None,
        scale: str | None = None,
        spikein_method: str | None = None,
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
        method:      Pileup method, e.g. ``"deeptools"``, ``"bamnado"``.
        scale:       Scaling method, e.g. ``"unscaled"``, ``"csaw"``, ``"spikein"``.
        spikein_method: Spike-in normalisation method, e.g. ``"orlando"``.
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
            df = df[df["method"] == method]
        if scale is not None:
            df = df[df["scale"] == scale]
        if spikein_method is not None:
            df = df[df["spikein_method"] == spikein_method]
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
        method: str | None = None,
        merged: bool | None = None,
        only_existing: bool = True,
    ) -> list[Path]:
        """Return peak BED paths matching the given filters."""
        df = self._peak_index
        if df.empty:
            return []

        targets = self._resolve_samples(sample, condition, antibody, group)
        if targets is not None:
            df = df[df["sample"].isin(targets)]
        if method is not None:
            df = df[df["method"] == method]
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
        method: str | None = None,
        scale: str | None = None,
        merged: bool | None = None,
        plot_type: str = "heatmap",
        only_existing: bool = True,
    ) -> list[Path]:
        """Return heatmap PDF paths.

        Parameters
        ----------
        plot_type: ``"heatmap"`` (default) or ``"metaplot"``.
        """
        heatmap_dir = self._output_dir / "heatmap"
        if not heatmap_dir.exists():
            return []

        paths = list(heatmap_dir.rglob(f"{plot_type}.pdf"))

        if method is not None:
            paths = [p for p in paths if f"/{method}/" in p.as_posix()]
        if scale is not None:
            paths = [p for p in paths if f"/{scale}/" in p.as_posix() or p.parent.name == scale]
        if merged is True:
            paths = [p for p in paths if "/merged/" in p.as_posix()]
        elif merged is False:
            paths = [p for p in paths if "/merged/" not in p.as_posix()]

        if only_existing:
            paths = [p for p in paths if p.exists()]
        return sorted(paths)

    def normalisation_factors(
        self, method: str | None = None
    ) -> list[Path]:
        """Return normalisation factor TSV files."""
        resources_dir = self._output_dir / "resources"
        if not resources_dir.exists():
            return []
        if method is not None:
            p = resources_dir / method / "normalisation_factors.tsv"
            return [p] if p.exists() else []
        return sorted(resources_dir.rglob("normalisation_factors.tsv"))

    def load_normalisation_factors(self, method: str | None = None) -> pd.DataFrame:
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

    def summary(self) -> pd.DataFrame:
        """Return a DataFrame summarising which file types are present."""
        rows = []
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
            rows.append({"File type": name, "Present": "✓" if present else "–"})
        return pd.DataFrame(rows)

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
