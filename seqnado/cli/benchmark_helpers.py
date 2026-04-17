"""Helpers for aggregating Snakemake benchmark TSV files into a report."""

from __future__ import annotations

import csv
import getpass
import hashlib
import math
import os
import re
from dataclasses import dataclass
from datetime import datetime
from html import escape
from pathlib import Path
from urllib.parse import quote

import pandas as pd

NUMERIC_COLUMNS = (
    "s",
    "h:m:s",
    "jobid",
    "max_rss",
    "max_vms",
    "max_uss",
    "max_pss",
    "io_in",
    "io_out",
    "mean_load",
    "cpu_time",
    "threads",
    "cpu_usage",
    "input_size_mb",
    "starttime",
    "endtime",
)


@dataclass(frozen=True)
class BenchmarkSummary:
    """Top-level metrics for a benchmark report."""

    jobs: int
    groups: int
    total_runtime_seconds: float
    total_cpu_time_seconds: float
    peak_rss_mb: float
    total_io_in: float
    total_io_out: float


TIMESTAMP_LINE_RE = re.compile(r"^\[(?P<timestamp>.+?)\]$")
JOB_START_RE = re.compile(r"^Job (?P<jobid>\d+): (?P<message>.+)$")
JOB_FINISH_RE = re.compile(r"^Finished jobid: (?P<jobid>\d+) \(Rule: (?P<rule>.+)\)$")
ALIGNMENT_READ_COUNT_RE = re.compile(
    r"^Aligned \((?P<rule>.+)\)\tReads mapped\t(?P<reads>\d+)$"
)
STEP_READ_COUNT_RE = re.compile(
    r"^Step \((?P<step>.+)\)\tReads Before\t(?P<before>\d+)\tReads After\t(?P<after>\d+)$"
)
STEP_READ_COUNT_LABEL_RE = re.compile(
    r"^(?P<step>[^\t]+)\tReads Before\t(?P<before>\d+)\tReads After\t(?P<after>\d+)$"
)
ENTITY_PATTERNS = (
    re.compile(r"\bsample (?P<entity>[^\s,;]+)"),
    re.compile(r"\bgroup (?P<entity>[^\s,;]+)"),
    re.compile(r"\bviewpoint (?P<entity>[^\s,;]+)"),
)
ASSAY_PREFIX_RE = re.compile(
    r"^(?P<assay>atac|chip|rna|meth|snp|cat|mcc|crispr|multiomics)_(?P<rule>.+)$"
)
ASSAY_COLORS = {
    "atac": "#c2410c",
    "chip": "#0f766e",
    "rna": "#2563eb",
    "meth": "#7c3aed",
    "snp": "#f6af2c",
    "cat": "#0891b2",
    "mcc": "#65a30d",
    "crispr": "#db2777",
    "multiomics": "#000000",
    "other": "#475569",
}
ASSAY_DISPLAY_NAMES = {
    "atac": "ATAC-seq",
    "chip": "ChIP-seq",
    "rna": "RNA-seq",
    "meth": "Methylation",
    "snp": "WGS (SNP)",
    "cat": "CUT&Tag",
    "mcc": "MCC",
    "crispr": "CRISPR",
    "multiomics": "Multiomics",
    "other": "Other",
}


def _assay_color(assay: str) -> str:
    """Return the configured color for an assay, falling back to `other`."""
    return ASSAY_COLORS.get(str(assay), ASSAY_COLORS["other"])


def _assay_display_name(assay: str) -> str:
    """Return the user-facing display name for an assay."""
    return ASSAY_DISPLAY_NAMES.get(str(assay), str(assay))


def _sample_color(sample: str) -> str:
    """Return a deterministic color for a sample series."""
    digest = hashlib.md5(sample.encode("utf-8")).hexdigest()
    hue = int(digest[:6], 16) % 360
    return f"hsl({hue}, 62%, 46%)"


def discover_benchmark_files(benchmark_dir: Path) -> list[Path]:
    """Return all benchmark TSV files under a benchmark directory."""
    return sorted(
        p
        for p in benchmark_dir.rglob("*.tsv")
        if p.is_file() and ".benchmark" in p.parts
    )


def discover_snakemake_logs(run_root: Path) -> list[Path]:
    """Return Snakemake log files under a workflow run root."""
    log_dir = run_root / ".snakemake" / "log"
    if not log_dir.exists():
        return []
    return sorted(p for p in log_dir.glob("*.snakemake.log") if p.is_file())


def discover_alignment_processing_logs(output_root: Path) -> list[Path]:
    """Return shared BAM-processing logs under a seqnado output root."""
    if not output_root.exists():
        return []
    return sorted(
        p for p in output_root.glob("*/qc/alignment_post_process/*.tsv") if p.is_file()
    )


def parse_snakemake_log_timeline(log_file: Path) -> pd.DataFrame:
    """Parse a Snakemake run log into a per-job timeline table."""
    rows: list[dict[str, object]] = []
    active_jobs: dict[str, dict[str, object]] = {}
    current_timestamp: datetime | None = None

    with open(log_file, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue

            timestamp_match = TIMESTAMP_LINE_RE.match(line)
            if timestamp_match:
                try:
                    current_timestamp = datetime.strptime(
                        timestamp_match.group("timestamp"), "%a %b %d %H:%M:%S %Y"
                    )
                except ValueError:
                    current_timestamp = None
                continue

            start_match = JOB_START_RE.match(line)
            if start_match and current_timestamp is not None:
                jobid = start_match.group("jobid")
                active_jobs[jobid] = {
                    "jobid": int(jobid),
                    "label": start_match.group("message"),
                    "entity": _extract_timeline_entity(start_match.group("message")),
                    "starttime": current_timestamp.timestamp(),
                }
                continue

            finish_match = JOB_FINISH_RE.match(line)
            if finish_match and current_timestamp is not None:
                jobid = finish_match.group("jobid")
                job = active_jobs.pop(jobid, None)
                if job is None:
                    job = {"jobid": int(jobid), "label": finish_match.group("rule")}
                rule_name = finish_match.group("rule")
                assay, display_rule = _split_rule_assay(rule_name)
                job["rule"] = rule_name
                job["assay"] = assay
                job["display_rule"] = display_rule
                job["endtime"] = current_timestamp.timestamp()
                rows.append(job)

    if not rows:
        return pd.DataFrame(
            columns=[
                "jobid",
                "rule",
                "display_rule",
                "assay",
                "label",
                "entity",
                "starttime",
                "endtime",
            ]
        )

    timeline = pd.DataFrame(rows)
    timeline = timeline[timeline["rule"] != "all"]
    return timeline.sort_values(["starttime", "endtime", "jobid"]).reset_index(
        drop=True
    )


def parse_snakemake_logs_timeline(log_files: list[Path]) -> pd.DataFrame:
    """Parse and merge multiple Snakemake logs, keeping the most recent rule/entity entry."""
    if not log_files:
        return pd.DataFrame(
            columns=[
                "jobid",
                "rule",
                "display_rule",
                "assay",
                "label",
                "entity",
                "starttime",
                "endtime",
            ]
        )

    timelines = [parse_snakemake_log_timeline(log_file) for log_file in log_files]
    timelines = [timeline for timeline in timelines if not timeline.empty]
    if not timelines:
        return pd.DataFrame(
            columns=[
                "jobid",
                "rule",
                "display_rule",
                "assay",
                "label",
                "entity",
                "starttime",
                "endtime",
            ]
        )

    merged = pd.concat(timelines, ignore_index=True)
    merged = merged.sort_values(["endtime", "starttime", "jobid"])
    merged = merged.drop_duplicates(subset=["rule", "entity"], keep="last")
    return merged.sort_values(["starttime", "endtime", "jobid"]).reset_index(drop=True)


def _extract_timeline_entity(message: str) -> str:
    """Extract a sample/group-like entity label from a Snakemake job message."""
    for pattern in ENTITY_PATTERNS:
        match = pattern.search(message)
        if match:
            return match.group("entity")
    return message


def _split_rule_assay(rule_name: str) -> tuple[str, str]:
    """Split an assay-prefixed rule name into assay and display rule parts."""
    match = ASSAY_PREFIX_RE.match(rule_name)
    if match:
        return match.group("assay"), match.group("rule")
    return "other", rule_name


def parse_alignment_processing_logs(output_root: Path) -> pd.DataFrame:
    """Parse shared BAM-processing logs into a read-count table."""
    rows: list[dict[str, object]] = []
    for log_file in discover_alignment_processing_logs(output_root):
        try:
            relative = log_file.relative_to(output_root)
        except ValueError:
            relative = log_file
        assay = relative.parts[0] if len(relative.parts) > 0 else "other"
        sample = log_file.stem
        try:
            tsv_df = pd.read_csv(log_file, sep="\t")
        except (OSError, pd.errors.EmptyDataError, pd.errors.ParserError):
            tsv_df = pd.DataFrame()

        if {"Step", "Reads Before", "Reads After"}.issubset(tsv_df.columns):
            for index, row in tsv_df.reset_index(drop=True).iterrows():
                step = str(row["Step"]).strip()
                if not step:
                    continue
                try:
                    read_count = int(float(row["Reads After"]))
                except (TypeError, ValueError):
                    continue
                sample_name = sample
                if "Sample" in tsv_df.columns:
                    raw_sample = str(row["Sample"]).strip()
                    if raw_sample:
                        sample_name = raw_sample
                rows.append(
                    {
                        "assay": assay,
                        "sample": sample_name,
                        "display_rule": step,
                        "read_count": read_count,
                        "sequence": index,
                    }
                )
            continue

        try:
            lines = log_file.read_text(encoding="utf-8").splitlines()
        except OSError:
            continue

        for index, raw_line in enumerate(lines):
            line = raw_line.strip()
            if not line:
                continue
            alignment_match = ALIGNMENT_READ_COUNT_RE.match(line)
            if alignment_match:
                rows.append(
                    {
                        "assay": assay,
                        "sample": sample,
                        "display_rule": f"aligned ({alignment_match.group('rule')})",
                        "read_count": int(alignment_match.group("reads")),
                        "sequence": index,
                    }
                )
                continue
            step_match = STEP_READ_COUNT_RE.match(line)
            if step_match:
                rows.append(
                    {
                        "assay": assay,
                        "sample": sample,
                        "display_rule": step_match.group("step"),
                        "read_count": int(step_match.group("after")),
                        "sequence": index,
                    }
                )
                continue
            step_label_match = STEP_READ_COUNT_LABEL_RE.match(line)
            if step_label_match:
                rows.append(
                    {
                        "assay": assay,
                        "sample": sample,
                        "display_rule": step_label_match.group("step"),
                        "read_count": int(step_label_match.group("after")),
                        "sequence": index,
                    }
                )

    if not rows:
        return pd.DataFrame(
            columns=["assay", "sample", "display_rule", "read_count", "sequence"]
        )

    return (
        pd.DataFrame(rows)
        .sort_values(["assay", "sample", "sequence"])
        .reset_index(drop=True)
    )


def _coerce_numeric_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Convert known benchmark metric columns to numeric where present."""
    for column in NUMERIC_COLUMNS:
        if column in df.columns:
            df[column] = pd.to_numeric(df[column], errors="coerce")
    return df


def _display_relative_benchmark_path(root_relative_path: Path) -> Path:
    """Return a cleaner path for reporting by stripping the .benchmark segment."""
    parts = list(root_relative_path.parts)
    if ".benchmark" not in parts:
        return root_relative_path

    benchmark_index = parts.index(".benchmark")
    cleaned_parts = parts[:benchmark_index] + parts[benchmark_index + 1 :]
    return Path(*cleaned_parts) if cleaned_parts else Path(root_relative_path.name)


def load_benchmark_table(benchmark_dir: Path) -> pd.DataFrame:
    """Load and annotate all benchmark TSV files from a directory tree."""
    rows: list[pd.DataFrame] = []

    for benchmark_file in discover_benchmark_files(benchmark_dir):
        relative_path = benchmark_file.relative_to(benchmark_dir)
        display_path = _display_relative_benchmark_path(relative_path)
        frame = pd.read_csv(benchmark_file, sep="\t")
        if frame.empty:
            continue

        frame = _coerce_numeric_columns(frame)
        frame["benchmark_file"] = str(display_path)
        frame["group"] = (
            str(display_path.parent) if display_path.parent != Path(".") else "root"
        )
        frame["job"] = benchmark_file.stem
        path_assay, path_display_rule = _split_benchmark_path(display_path)
        if "rule_name" in frame.columns:
            parsed_rule_names = (
                frame["rule_name"].fillna("").astype(str).map(_split_rule_assay)
            )
            frame["assay"] = parsed_rule_names.map(
                lambda parts: parts[0] if parts[0] != "other" else path_assay
            )
            frame["display_rule"] = parsed_rule_names.map(
                lambda parts: (
                    parts[1] if parts[1] and parts[1] != "root" else path_display_rule
                )
            )
        else:
            frame["assay"] = path_assay
            frame["display_rule"] = path_display_rule
        rows.append(frame)

    if not rows:
        return pd.DataFrame()

    combined = pd.concat(rows, ignore_index=True)
    if "s" in combined.columns:
        combined = combined.sort_values(
            ["s", "benchmark_file"], ascending=[False, True]
        )
    return combined.reset_index(drop=True)


def _split_benchmark_path(display_path: Path) -> tuple[str, str]:
    """Infer assay and rule-like grouping from a benchmark path."""
    parts = list(display_path.parts)
    if not parts:
        return "other", "root"

    assay = parts[0] if parts[0] in ASSAY_COLORS else "other"
    rule_parts = parts[1:-1] if assay != "other" else parts[:-1]
    display_rule = "/".join(rule_parts) if rule_parts else "root"
    return assay, display_rule


def _format_rule_label(rule: str) -> str:
    """Convert a rule/path identifier into a human-readable label."""
    if not rule or rule == "root":
        return "root"
    return rule.replace("_", " ").replace("/", " / ")


def _shortest_unique_rule_labels(rules: list[str]) -> dict[str, str]:
    """Build the shortest unique human-readable suffix label for each rule path."""
    unique_rules = list(dict.fromkeys(rules))
    labels: dict[str, str] = {}

    for rule in unique_rules:
        parts = rule.split("/") if rule else ["root"]
        chosen = rule
        for suffix_len in range(1, len(parts) + 1):
            candidate = "/".join(parts[-suffix_len:])
            collisions = 0
            for other in unique_rules:
                other_parts = other.split("/") if other else ["root"]
                if (
                    len(other_parts) >= suffix_len
                    and "/".join(other_parts[-suffix_len:]) == candidate
                ):
                    collisions += 1
            if collisions == 1:
                chosen = candidate
                break
        labels[rule] = _format_rule_label(chosen)

    return labels


def summarize_benchmarks(df: pd.DataFrame) -> BenchmarkSummary:
    """Compute top-level summary metrics for a benchmark table."""
    if df.empty:
        return BenchmarkSummary(0, 0, 0.0, 0.0, 0.0, 0.0, 0.0)

    return BenchmarkSummary(
        jobs=len(df),
        groups=df["group"].nunique() if "group" in df.columns else 0,
        total_runtime_seconds=float(df["s"].fillna(0).sum())
        if "s" in df.columns
        else 0.0,
        total_cpu_time_seconds=float(df["cpu_time"].fillna(0).sum())
        if "cpu_time" in df.columns
        else 0.0,
        peak_rss_mb=float(df["max_rss"].fillna(0).max())
        if "max_rss" in df.columns
        else 0.0,
        total_io_in=float(df["io_in"].fillna(0).sum())
        if "io_in" in df.columns
        else 0.0,
        total_io_out=float(df["io_out"].fillna(0).sum())
        if "io_out" in df.columns
        else 0.0,
    )


def format_seconds(seconds: float) -> str:
    """Format seconds in a compact human-readable form."""
    total_seconds = int(round(seconds))
    hours, remainder = divmod(total_seconds, 3600)
    minutes, secs = divmod(remainder, 60)
    return f"{hours:d}:{minutes:02d}:{secs:02d}"


def _format_metric(value: float | int | None, decimals: int = 2) -> str:
    """Render a numeric metric, preserving blanks for missing values."""
    if value is None or pd.isna(value):
        return ""
    return f"{float(value):,.{decimals}f}"


def format_compact_number(value: float | int | None, decimals: int = 1) -> str:
    """Render large numeric values using compact suffixes such as k, M, and B."""
    if value is None or pd.isna(value):
        return ""

    magnitude = abs(float(value))
    suffixes = (
        (1_000_000_000_000, "T"),
        (1_000_000_000, "B"),
        (1_000_000, "M"),
        (1_000, "k"),
    )
    for threshold, suffix in suffixes:
        if magnitude >= threshold:
            suffix_decimals = 1 if suffix == "k" else 2
            return f"{float(value) / threshold:.{suffix_decimals}f}{suffix}"
    return _format_metric(value, decimals=0 if float(value).is_integer() else decimals)


def format_bytes(value: float | int | None) -> str:
    """Render byte counts in human-readable binary units."""
    if value is None or pd.isna(value):
        return ""

    size = float(value)
    units = ("B", "KiB", "MiB", "GiB", "TiB")
    unit_index = 0
    while abs(size) >= 1024 and unit_index < len(units) - 1:
        size /= 1024
        unit_index += 1

    decimals = 0 if unit_index == 0 else (2 if abs(size) < 10 else 1)
    return f"{size:.{decimals}f} {units[unit_index]}"


def _csv_download_link(df: pd.DataFrame, filename: str) -> str:
    """Render a CSV download link for a dataframe."""
    if df.empty:
        return ""

    csv_text = df.to_csv(index=False)
    href = "data:text/csv;charset=utf-8," + quote(csv_text)
    return f"<a class='export-link' download='{escape(filename)}' href='{href}'>Export CSV</a>"


def _plot_description_html(description: str) -> str:
    """Render a short description beneath a plot export link."""
    return f"<p class='plot-description'>{escape(description)}</p>"


def _expanded_metric_description(value_column: str, value_label: str) -> str:
    """Return a human-readable metric description with acronyms expanded."""
    descriptions = {
        "s": "runtime in seconds",
        "max_rss": "maximum resident set size (peak memory use) in megabytes",
        "max_vms": "maximum virtual memory size in megabytes",
        "max_uss": "maximum unique set size in megabytes",
        "max_pss": "maximum proportional set size in megabytes",
        "io_in": "input/output read volume",
        "io_out": "input/output write volume",
        "io_signed": "input/output volume with reads shown as negative values and writes shown as positive values",
        "cpu_usage": "CPU usage percentage",
        "input_size_mb": "input size in megabytes",
        "read_count": "read counts extracted from shared BAM processing logs",
        "output_size_bytes": "output size in bytes",
    }
    return descriptions.get(value_column, value_label.lower())


def _metric_formatter(value_column: str):
    """Return the formatter function appropriate for a metric column."""
    if value_column in {"io_in", "io_out"}:
        return format_compact_number
    if value_column == "output_size_bytes":
        return format_bytes
    return _format_metric


def _metric_available(df: pd.DataFrame, column: str) -> bool:
    """Return whether a metric column is present with at least one non-null value."""
    return column in df.columns and df[column].notna().any()


def _collapsible_section_html(
    section_id: str, title: str, content: str, open_by_default: bool = False
) -> str:
    """Wrap report content in a collapsible section with a stable anchor."""
    open_attr = " open" if open_by_default else ""
    return (
        f"<details class='report-section' id='{escape(section_id)}'{open_attr}>"
        f"<summary><span>{escape(title)}</span></summary>"
        f"<div class='section section-body'>{content}</div>"
        "</details>"
    )


def _assay_swatch_html(assay: str) -> str:
    """Render a colored assay swatch."""
    return f"<span class='legend-swatch' style='background:{escape(_assay_color(assay))}'></span>"


def _report_run_root(benchmark_dir: Path) -> Path:
    """Return the workflow run root for a benchmark directory."""
    return (
        benchmark_dir.parent
        if benchmark_dir.name == "seqnado_output"
        else benchmark_dir
    )


def _assay_chip_legend_html(assays: list[str]) -> str:
    """Render a compact non-interactive legend for assay colors."""
    if not assays:
        return ""
    items = "".join(
        (
            "<span class='legend-item'>"
            f"{_assay_swatch_html(assay)}"
            f"<span class='legend-label'>{escape(_assay_display_name(assay))}</span>"
            "</span>"
        )
        for assay in assays
    )
    return f"<div class='legend chart-legend'>{items}</div>"


def _compute_assay_sample_counts(
    benchmark_dir: Path, timeline_df: pd.DataFrame
) -> dict[str, int]:
    """Compute per-assay sample counts, preferring metadata CSVs over inferred timeline entities."""
    run_root = _report_run_root(benchmark_dir)
    metadata_paths = sorted(run_root.glob("metadata_*.csv"))
    counts: dict[str, int] = {}

    for metadata_path in metadata_paths:
        try:
            with open(metadata_path, newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
        except OSError:
            continue
        if not rows:
            continue
        assay_name = str(rows[0].get("assay", "")).strip().lower()
        if not assay_name:
            assay_name = metadata_path.stem.removeprefix("metadata_").lower()
        sample_ids = {
            str(row.get("sample_id", "")).strip()
            for row in rows
            if str(row.get("sample_id", "")).strip()
        }
        if sample_ids:
            counts[assay_name] = len(sample_ids)

    if counts:
        return counts

    if timeline_df.empty or not {"assay", "entity"}.issubset(timeline_df.columns):
        return {}

    counts_df = timeline_df.loc[:, ["assay", "entity"]].dropna().copy()
    counts_df["entity"] = counts_df["entity"].astype(str).str.strip()
    counts_df = counts_df[counts_df["entity"] != ""]
    if counts_df.empty:
        return {}

    return (
        counts_df.groupby("assay")["entity"]
        .nunique()
        .sort_values(ascending=False)
        .to_dict()
    )


def _overview_assay_counts_html(benchmark_dir: Path, timeline_df: pd.DataFrame) -> str:
    """Render sample counts per assay for the overview section."""
    counts = _compute_assay_sample_counts(benchmark_dir, timeline_df)
    if not counts:
        return ""

    items = "".join(
        (
            "<div class='assay-count-chip'>"
            f"{_assay_swatch_html(assay)}"
            f"<span class='assay-count-label'>{escape(_assay_display_name(assay))}</span>"
            f"<span class='assay-count-value'>{count}</span>"
            "</div>"
        )
        for assay, count in counts.items()
    )
    return (
        "<div class='overview-assay-counts'>"
        "<div class='overview-subtitle'>Samples Per Assay</div>"
        f"<div class='assay-count-grid'>{items}</div>"
        "</div>"
    )


def _seqnado_version() -> str:
    """Return the local SeqNado version string when available."""
    version_file = Path(__file__).resolve().parent.parent / "_version.py"
    if not version_file.exists():
        return "unknown"
    match = re.search(
        r"__version__ = version = '([^']+)'", version_file.read_text(encoding="utf-8")
    )
    return match.group(1) if match else "unknown"


def _overview_provenance_html() -> str:
    """Render report provenance details for the overview section."""
    generated_at = datetime.now().astimezone().strftime("%Y-%m-%d %H:%M:%S %Z")
    user_name = (
        os.environ.get("USER")
        or os.environ.get("USERNAME")
        or os.environ.get("LOGNAME")
        or getpass.getuser()
        or "unknown"
    )
    return (
        "<div class='overview-provenance'>"
        "<div class='overview-subtitle'>Report Generation</div>"
        "<div class='provenance-line'>"
        f"<span>User <strong>{escape(user_name)}</strong> generated this report with SeqNado version <strong>{escape(_seqnado_version())}</strong> on <strong>{escape(generated_at)}</strong>.</span>"
        "</div>"
        "</div>"
    )


def _table_html(df: pd.DataFrame, columns: list[str], title: str) -> str:
    """Render a dataframe slice as a small HTML table."""
    safe = df.loc[:, [c for c in columns if c in df.columns]].copy()
    if safe.empty:
        return f"<h2>{escape(title)}</h2><p>No data available.</p>"

    display = safe.rename(
        columns={
            "benchmark_file": "Benchmark File",
            "group": "Group",
            "job": "Job",
            "s": "Runtime (s)",
            "cpu_time": "CPU Time (s)",
            "max_rss": "Max RSS (MB)",
            "io_in": "I/O In",
            "io_out": "I/O Out",
            "mean_load": "Mean Load",
        }
    )
    for column in ("Runtime (s)", "CPU Time (s)", "Max RSS (MB)", "Mean Load"):
        if column in display.columns:
            display[column] = display[column].map(_format_metric)
    for column in ("I/O In", "I/O Out"):
        if column in display.columns:
            display[column] = display[column].map(format_compact_number)
    if "Output Size" in display.columns:
        display["Output Size"] = display["Output Size"].map(format_bytes)

    return f"<h2>{escape(title)}</h2>{display.to_html(index=False, classes='report-table', border=0, escape=True)}"


def _bar_plot_html(
    df: pd.DataFrame,
    label_column: str,
    value_column: str,
    title: str,
    value_label: str,
    top_n: int = 10,
) -> str:
    """Render a lightweight horizontal bar chart as inline SVG."""
    if df.empty or label_column not in df.columns or value_column not in df.columns:
        return f"<h2>{escape(title)}</h2><p>No data available.</p>"

    plot_df = (
        df.loc[:, [label_column, value_column]]
        .dropna(subset=[value_column])
        .sort_values(value_column, ascending=False)
        .head(top_n)
    )
    if plot_df.empty:
        return f"<h2>{escape(title)}</h2><p>No data available.</p>"

    max_value = float(plot_df[value_column].max())
    if max_value <= 0:
        return f"<h2>{escape(title)}</h2><p>No positive values available.</p>"

    width = 900
    left_margin = 260
    right_margin = 110
    top_margin = 28
    row_height = 28
    bar_height = 18
    bottom_margin = 46
    height = top_margin + len(plot_df) * row_height + bottom_margin
    chart_width = width - left_margin - right_margin
    tick_count = 5

    svg_parts = [
        f"<h2>{escape(title)}</h2>",
        _csv_download_link(plot_df, f"{title.lower().replace(' ', '_')}.csv"),
        _plot_description_html(
            "Total output size per assay, measured from files in each assay directory under seqnado_output and excluding .benchmark files."
            if value_column == "output_size_bytes" and label_column == "assay"
            else f"Top {min(len(plot_df), top_n)} entries ranked by {value_label.lower()}."
        ),
        f"<svg viewBox='0 0 {width} {height}' class='plot' role='img' aria-label='{escape(title)}'>",
    ]

    for idx, row in enumerate(plot_df.itertuples(index=False), start=0):
        label = str(getattr(row, label_column))
        label_text = _assay_display_name(label) if label_column == "assay" else label
        value = float(getattr(row, value_column))
        y = top_margin + idx * row_height
        bar_width = 0 if max_value == 0 else chart_width * (value / max_value)
        bar_fill = _assay_color(str(label)) if label_column == "assay" else "#b45309"
        display_value = _metric_formatter(value_column)(value)
        svg_parts.extend(
            [
                f"<text x='{left_margin - 10}' y='{y + 13}' text-anchor='end' class='plot-label'>{escape(label_text)}</text>",
                f"<rect x='{left_margin}' y='{y}' width='{chart_width}' height='{bar_height}' rx='6' ry='6' class='plot-track'></rect>",
                f"<rect x='{left_margin}' y='{y}' width='{bar_width:.2f}' height='{bar_height}' rx='6' ry='6' class='plot-bar' style='fill:{bar_fill}'><title>{escape(label_text)} | {escape(display_value)}</title></rect>",
                f"<text x='{left_margin + chart_width + 8}' y='{y + 13}' class='plot-value'>{escape(display_value)}</text>",
            ]
        )

    axis_y = height - 24
    svg_parts.append(
        f"<line x1='{left_margin}' y1='{axis_y}' x2='{left_margin + chart_width}' y2='{axis_y}' class='axis-line'></line>"
    )
    for tick_index in range(tick_count + 1):
        fraction = tick_index / tick_count
        tick_x = left_margin + chart_width * fraction
        tick_value = max_value * fraction
        tick_label = _metric_formatter(value_column)(tick_value)
        svg_parts.extend(
            [
                f"<line x1='{tick_x:.2f}' y1='{axis_y}' x2='{tick_x:.2f}' y2='{axis_y + 6}' class='axis-line'></line>",
                f"<text x='{tick_x:.2f}' y='{axis_y + 18}' text-anchor='middle' class='plot-axis'>{escape(tick_label)}</text>",
            ]
        )
    svg_parts.append(
        f"<text x='{left_margin}' y='{height - 34}' class='plot-axis'>{escape(value_label)}</text>"
    )
    svg_parts.append("</svg>")
    return "".join(svg_parts)


def _box_plot_html(
    df: pd.DataFrame,
    value_column: str,
    title: str,
    value_label: str,
    top_n: int = 12,
) -> str:
    """Render an assay-colored horizontal box plot grouped by derived rule."""
    required = {"assay", "display_rule", value_column}
    if df.empty or not required.issubset(df.columns):
        return f"<h2>{escape(title)}</h2><p>No data available.</p>"

    plot_df = df.loc[:, ["assay", "display_rule", value_column]].copy()
    plot_df[value_column] = pd.to_numeric(plot_df[value_column], errors="coerce")
    plot_df = plot_df.dropna(subset=[value_column])
    if plot_df.empty:
        return f"<h2>{escape(title)}</h2><p>No data available.</p>"

    group_summary = (
        plot_df.groupby(["assay", "display_rule"], dropna=False)[value_column]
        .agg(["count", "median"])
        .reset_index()
        .sort_values(["median", "count"], ascending=[False, False])
        .head(top_n)
    )
    if group_summary.empty:
        return f"<h2>{escape(title)}</h2><p>No data available.</p>"

    selected = plot_df.merge(
        group_summary[["assay", "display_rule"]],
        on=["assay", "display_rule"],
        how="inner",
    )
    order = list(group_summary.itertuples(index=False, name=None))
    stats_rows: list[dict[str, object]] = []
    for assay, display_rule, _, _ in order:
        values = selected.loc[
            (selected["assay"] == assay) & (selected["display_rule"] == display_rule),
            value_column,
        ].sort_values()
        q1 = float(values.quantile(0.25))
        median = float(values.quantile(0.5))
        q3 = float(values.quantile(0.75))
        min_v = float(values.min())
        max_v = float(values.max())
        stats_rows.append(
            {
                "assay": assay,
                "display_rule": display_rule,
                "label": f"{assay} | {display_rule}",
                "count": int(values.count()),
                "min": min_v,
                "q1": q1,
                "median": median,
                "q3": q3,
                "max": max_v,
            }
        )

    stats_df = pd.DataFrame(stats_rows)
    rule_order = (
        stats_df.groupby("display_rule", sort=False)["median"]
        .max()
        .sort_values(ascending=False)
        .index.tolist()
    )
    rule_positions = {rule: idx for idx, rule in enumerate(rule_order)}
    rule_labels = _shortest_unique_rule_labels(rule_order)
    present_assays = [
        assay
        for assay in ASSAY_COLORS
        if assay != "other" and assay in stats_df["assay"].unique()
    ]
    if "other" in stats_df["assay"].unique():
        present_assays.append("other")
    assay_spacing = 14
    assay_order = {assay: idx for idx, assay in enumerate(present_assays)}
    max_assays_per_rule = int(stats_df.groupby("display_rule")["assay"].nunique().max())
    max_value = float(stats_df["max"].max())
    min_value = float(stats_df["min"].min())
    center_zero = value_column == "io_signed"
    if center_zero:
        max_abs_value = max(abs(min_value), abs(max_value), 1.0)
        use_log_scale = False
    else:
        if max_value <= 0:
            max_value = 1.0
        positive_mins = stats_df.loc[stats_df["min"] > 0, "min"]
        min_positive = float(positive_mins.min()) if not positive_mins.empty else 0.0
        use_log_scale = bool(min_positive > 0 and max_value / min_positive >= 100)

    width = 1040
    left_margin = 300
    right_margin = 70
    top_margin = 28
    row_height = max(52, 24 + (max_assays_per_rule + 1) * assay_spacing)
    box_height = 12
    bottom_margin = 52
    chart_width = width - left_margin - right_margin
    height = top_margin + len(rule_order) * row_height + bottom_margin
    tick_count = 5

    def scale_x(value: float) -> float:
        if center_zero:
            return left_margin + chart_width * (
                (value + max_abs_value) / (2 * max_abs_value)
            )
        if use_log_scale:
            clipped = max(value, min_positive)
            domain_min = math.log10(min_positive)
            domain_max = math.log10(max_value)
            if domain_max == domain_min:
                return left_margin
            return left_margin + chart_width * (
                (math.log10(clipped) - domain_min) / (domain_max - domain_min)
            )
        return left_margin + chart_width * (value / max_value)

    value_formatter = _metric_formatter(value_column)

    svg_parts = [
        f"<h2>{escape(title)}</h2>",
        _csv_download_link(selected, f"{title.lower().replace(' ', '_')}.csv"),
        _plot_description_html(
            f"Distribution of {_expanded_metric_description(value_column, value_label)} for the highest-signal rule and assay combinations. Each box summarizes all matching benchmark records with min, lower quartile, median, upper quartile, and max."
            + (
                " Negative values indicate reads and positive values indicate writes."
                if center_zero
                else ""
            )
            + (
                " The x-axis uses a logarithmic scale because the value range spans multiple orders of magnitude."
                if use_log_scale
                else ""
            )
        ),
        _assay_chip_legend_html(present_assays),
        f"<svg viewBox='0 0 {width} {height}' class='plot' role='img' aria-label='{escape(title)}'>",
    ]

    if center_zero:
        tick_values = [max_abs_value * fraction for fraction in (-1, -0.5, 0, 0.5, 1)]
    elif use_log_scale:
        start_exp = math.floor(math.log10(min_positive))
        end_exp = math.ceil(math.log10(max_value))
        tick_values = [10**exp for exp in range(start_exp, end_exp + 1)]
        if min_positive not in tick_values:
            tick_values = [min_positive] + tick_values
        if max_value not in tick_values:
            tick_values.append(max_value)
        tick_values = sorted({v for v in tick_values if min_positive <= v <= max_value})
    else:
        tick_values = [
            max_value * (tick_index / tick_count)
            for tick_index in range(tick_count + 1)
        ]

    for tick_value in tick_values:
        tick_x = scale_x(float(tick_value))
        svg_parts.append(
            f"<line x1='{tick_x:.2f}' y1='{top_margin - 4}' x2='{tick_x:.2f}' y2='{height - bottom_margin + 8}' class='plot-grid-line'></line>"
        )

    for rule, row_index in rule_positions.items():
        band_y = top_margin + row_index * row_height
        line_y = band_y + row_height - 10
        svg_parts.extend(
            [
                f"<rect x='{left_margin}' y='{band_y:.2f}' width='{chart_width}' height='{row_height:.2f}' class='plot-row-band'></rect>",
                f"<text x='{left_margin - 10}' y='{line_y + 4:.2f}' text-anchor='end' class='plot-label'>{escape(rule_labels.get(str(rule), str(rule)))}</text>",
                f"<line x1='{left_margin}' y1='{line_y:.2f}' x2='{left_margin + chart_width}' y2='{line_y:.2f}' class='plot-track-line'></line>",
            ]
        )

    for idx, row in enumerate(stats_df.itertuples(index=False), start=0):
        band_y = top_margin + rule_positions[str(row.display_rule)] * row_height
        line_y = band_y + row_height - 10
        assay_index = assay_order.get(str(row.assay), 0)
        y = line_y - box_height - ((assay_index + 1) * assay_spacing)
        color = _assay_color(str(row.assay))
        min_x = scale_x(float(row.min))
        q1_x = scale_x(float(row.q1))
        median_x = scale_x(float(row.median))
        q3_x = scale_x(float(row.q3))
        max_x = scale_x(float(row.max))
        label = f"{row.display_rule} | {_assay_display_name(str(row.assay))}"
        tooltip = (
            f"{label} | n={row.count} | min={value_formatter(row.min)} | "
            f"q1={value_formatter(row.q1)} | median={value_formatter(row.median)} | "
            f"q3={value_formatter(row.q3)} | max={value_formatter(row.max)}"
        )
        svg_parts.extend(
            [
                f"<line x1='{min_x:.2f}' y1='{y + box_height / 2:.2f}' x2='{max_x:.2f}' y2='{y + box_height / 2:.2f}' class='box-whisker'></line>",
                f"<line x1='{min_x:.2f}' y1='{y + 2:.2f}' x2='{min_x:.2f}' y2='{y + box_height - 2:.2f}' class='box-whisker-cap'></line>",
                f"<line x1='{max_x:.2f}' y1='{y + 2:.2f}' x2='{max_x:.2f}' y2='{y + box_height - 2:.2f}' class='box-whisker-cap'></line>",
                f"<rect x='{q1_x:.2f}' y='{y:.2f}' width='{max(q3_x - q1_x, 3):.2f}' height='{box_height}' rx='6' ry='6' fill='{color}' class='box-rect'><title>{escape(tooltip)}</title></rect>",
                f"<line x1='{median_x:.2f}' y1='{y:.2f}' x2='{median_x:.2f}' y2='{y + box_height:.2f}' class='box-median'></line>",
                f"<circle cx='{median_x:.2f}' cy='{y + box_height / 2:.2f}' r='2.4' class='box-median-dot'></circle>",
            ]
        )

    axis_y = height - 24
    svg_parts.append(
        f"<line x1='{left_margin}' y1='{axis_y}' x2='{left_margin + chart_width}' y2='{axis_y}' class='axis-line'></line>"
    )
    for tick_value in tick_values:
        tick_x = scale_x(float(tick_value))
        tick_label = value_formatter(tick_value)
        svg_parts.extend(
            [
                f"<line x1='{tick_x:.2f}' y1='{axis_y}' x2='{tick_x:.2f}' y2='{axis_y + 6}' class='axis-line'></line>",
                f"<text x='{tick_x:.2f}' y='{axis_y + 18}' text-anchor='middle' class='plot-axis'>{escape(tick_label)}</text>",
            ]
        )
    if center_zero:
        zero_x = scale_x(0)
        svg_parts.append(
            f"<line x1='{zero_x:.2f}' y1='{top_margin - 4}' x2='{zero_x:.2f}' y2='{height - bottom_margin + 8}' class='plot-zero-line'></line>"
        )
    svg_parts.append(
        f"<text x='{left_margin}' y='{height - 34}' class='plot-axis'>{escape(value_label + (' (log scale)' if use_log_scale else ''))}</text>"
    )
    svg_parts.append("</svg>")
    return "".join(svg_parts)


def _read_count_line_plot_html(read_counts_df: pd.DataFrame) -> str:
    """Render a per-sample read-count trajectory plot."""
    required = {"sample", "display_rule", "read_count", "sequence"}
    if read_counts_df.empty or not required.issubset(read_counts_df.columns):
        return "<h2>Read Counts</h2><p>No data available.</p>"

    plot_df = read_counts_df.loc[
        :, ["sample", "display_rule", "read_count", "sequence"]
    ].copy()
    plot_df["read_count"] = pd.to_numeric(plot_df["read_count"], errors="coerce")
    plot_df["sequence"] = pd.to_numeric(plot_df["sequence"], errors="coerce")
    plot_df = plot_df.dropna(subset=["read_count", "sequence"])
    if plot_df.empty:
        return "<h2>Read Counts</h2><p>No data available.</p>"

    plot_df = plot_df.sort_values(["sample", "sequence", "display_rule"]).reset_index(
        drop=True
    )
    plot_df["step_occurrence"] = (
        plot_df.groupby(["sample", "display_rule"]).cumcount() + 1
    )

    canonical_intermediate_order = {
        ("Rename Aligned BAM", 1): 0,
        ("QNAME Sort", 1): 2,
        ("Spike-in Filter", 1): 3,
        ("Spike-in Split", 1): 5,
        ("Move Reference BAM", 1): 6,
        ("QNAME Sort", 2): 8,
        ("Blacklist", 1): 10,
        ("Remove Duplicates", 1): 11,
        ("ATAC Shift", 1): 12,
        ("Filter", 1): 13,
    }
    aligned_mask = (
        plot_df["display_rule"].astype(str).str.lower().str.startswith("aligned (")
    )
    finalise_mask = (
        plot_df["display_rule"].astype(str).str.strip().str.lower() == "finalise"
    )

    intermediate_df = plot_df.loc[
        ~aligned_mask & ~finalise_mask, ["display_rule", "sequence", "step_occurrence"]
    ].copy()
    step_occurrence_max = (
        intermediate_df.groupby("display_rule", dropna=False)["step_occurrence"]
        .max()
        .to_dict()
        if not intermediate_df.empty
        else {}
    )
    observed_intermediate_order = (
        intermediate_df.groupby("display_rule", dropna=False)["sequence"]
        .median()
        .sort_values()
        .index.tolist()
        if not intermediate_df.empty
        else []
    )
    axis_slots: list[tuple[str, int, str]] = [("__aligned__", 1, "Aligned")]
    observed_slots: list[tuple[str, int, str]] = []
    for step in observed_intermediate_order:
        occurrences = int(step_occurrence_max.get(step, 1))
        for occurrence in range(1, occurrences + 1):
            label = f"{step} {occurrence}" if occurrences > 1 else str(step)
            observed_slots.append((str(step), occurrence, label))
    observed_slots = sorted(
        observed_slots,
        key=lambda slot: (
            canonical_intermediate_order.get(
                (slot[0], slot[1]), len(canonical_intermediate_order)
            ),
            observed_intermediate_order.index(slot[0]),
            slot[1],
        ),
    )
    axis_slots.extend(observed_slots)
    axis_slots.append(("Finalise", 1, "Finalise"))

    step_positions = {
        (step, occurrence): idx
        for idx, (step, occurrence, _label) in enumerate(axis_slots)
    }

    def _step_position(row: pd.Series) -> int:
        step = str(row["display_rule"])
        if step.lower().startswith("aligned ("):
            return step_positions[("__aligned__", 1)]
        if step.strip().lower() == "finalise":
            return step_positions[("Finalise", 1)]
        return step_positions.get(
            (step, int(row["step_occurrence"])), step_positions[("Finalise", 1)] - 1
        )

    plot_df["step_position"] = plot_df.apply(_step_position, axis=1)
    plot_df["reads_millions"] = plot_df["read_count"] / 1_000_000
    sample_order = sorted(
        plot_df["sample"].dropna().astype(str).unique().tolist(), key=str.casefold
    )

    unique_positions = sorted(set(step_positions.values()))
    width = max(960, 180 + max(len(unique_positions), 1) * 110)
    height = 520
    left_margin = 80
    right_margin = 28
    top_margin = 24
    bottom_margin = 120
    chart_width = width - left_margin - right_margin
    chart_height = height - top_margin - bottom_margin
    max_reads_m = float(plot_df["reads_millions"].max()) if not plot_df.empty else 1.0
    if max_reads_m <= 0:
        max_reads_m = 1.0
    tick_count = 5

    def scale_x(position: int) -> float:
        if len(unique_positions) <= 1:
            return left_margin + chart_width / 2
        max_position = max(unique_positions)
        if max_position <= 0:
            return left_margin + chart_width / 2
        return left_margin + chart_width * (position / max_position)

    def scale_y(value: float) -> float:
        return top_margin + chart_height - (chart_height * (value / max_reads_m))

    legend_items = "".join(
        (
            f"<button type='button' class='chip legend-toggle is-active' data-sample='{escape(str(sample))}'>"
            f"<span class='chip-swatch' style='background:{_sample_color(str(sample))}'></span>"
            f"{escape(str(sample))}"
            "</button>"
        )
        for sample in sample_order
    )

    svg_parts = [
        "<h2>Read Counts</h2>",
        _csv_download_link(plot_df, "read_count_trajectory.csv"),
        _plot_description_html(
            "Read-count trajectories across BAM processing steps. Reads are measured in millions are per fragment for paired-end data and per read for single-end data. Hover over points to see exact read counts and step details. Toggle samples on/off by clicking their legend chips."
        ),
        f"<div class='legend read-count-legend'>{legend_items}</div>",
        "<div class='plot-tooltip' id='read-count-tooltip' hidden></div>",
        f"<svg viewBox='0 0 {width} {height}' class='plot' id='read-count-chart' role='img' aria-label='Read Count Trajectory'>",
    ]

    for tick_index in range(tick_count + 1):
        tick_value = max_reads_m * (tick_index / tick_count)
        tick_y = scale_y(tick_value)
        svg_parts.extend(
            [
                f"<line x1='{left_margin}' y1='{tick_y:.2f}' x2='{left_margin + chart_width}' y2='{tick_y:.2f}' class='plot-grid-line'></line>",
                f"<text x='{left_margin - 10}' y='{tick_y + 4:.2f}' text-anchor='end' class='plot-axis'>{tick_value:,.2f}</text>",
            ]
        )

    position_labels = {
        idx: label for idx, (_step, _occurrence, label) in enumerate(axis_slots)
    }

    for position in unique_positions:
        tick_x = scale_x(position)
        label = escape(position_labels.get(position, ""))
        svg_parts.extend(
            [
                f"<line x1='{tick_x:.2f}' y1='{top_margin}' x2='{tick_x:.2f}' y2='{top_margin + chart_height}' class='plot-grid-line'></line>",
                f"<text x='{tick_x:.2f}' y='{height - 62}' text-anchor='end' transform='rotate(-35 {tick_x:.2f} {height - 62})' class='plot-axis'>{label}</text>",
            ]
        )

    svg_parts.append(
        f"<line x1='{left_margin}' y1='{top_margin + chart_height}' x2='{left_margin + chart_width}' y2='{top_margin + chart_height}' class='axis-line'></line>"
    )
    svg_parts.append(
        f"<line x1='{left_margin}' y1='{top_margin}' x2='{left_margin}' y2='{top_margin + chart_height}' class='axis-line'></line>"
    )
    svg_parts.append(
        f"<text x='{left_margin}' y='{height - 18}' class='plot-axis'>Workflow Step</text>"
    )
    svg_parts.append(
        f"<text x='18' y='{top_margin + chart_height / 2:.2f}' transform='rotate(-90 18 {top_margin + chart_height / 2:.2f})' class='plot-axis'>Reads (millions)</text>"
    )

    for sample in sample_order:
        sample_df = plot_df[plot_df["sample"] == sample].sort_values(
            ["step_position", "sequence"]
        )
        if sample_df.empty:
            continue
        points = " ".join(
            f"{scale_x(int(row.step_position)):.2f},{scale_y(float(row.reads_millions)):.2f}"
            for row in sample_df.itertuples(index=False)
        )
        color = _sample_color(str(sample))
        svg_parts.append(
            f"<polyline class='read-count-series' data-sample='{escape(str(sample))}' fill='none' stroke='{color}' stroke-width='2.5' stroke-linecap='round' stroke-linejoin='round' points='{points}'><title>{escape(str(sample))}</title></polyline>"
        )
        for row in sample_df.itertuples(index=False):
            point_x = scale_x(int(row.step_position))
            point_y = scale_y(float(row.reads_millions))
            svg_parts.append(
                f"<circle class='read-count-point' data-sample='{escape(str(sample))}' data-step='{escape(str(row.display_rule))}' data-reads='{int(row.read_count):,}' data-reads-millions='{row.reads_millions:,.3f}' cx='{point_x:.2f}' cy='{point_y:.2f}' r='3.5' fill='{color}'><title>{escape(str(sample))} | {escape(str(row.display_rule))} | {row.reads_millions:,.2f}M reads</title></circle>"
            )

    svg_parts.append("</svg>")
    return "".join(svg_parts)


def _gantt_plot_html(df: pd.DataFrame, title: str) -> str:
    """Render a Gantt-style timeline with rule rows and sample/group bars."""
    required = {"display_rule", "assay", "entity", "starttime", "endtime"}
    if df.empty or not required.issubset(df.columns):
        return f"<h2>{escape(title)}</h2><p>Gantt chart unavailable: no Snakemake run timeline could be parsed.</p>"

    plot_df = df.loc[
        :, ["display_rule", "assay", "entity", "starttime", "endtime"]
    ].dropna()
    if plot_df.empty:
        return f"<h2>{escape(title)}</h2><p>Gantt chart unavailable: start/end timestamps are missing.</p>"

    plot_df = plot_df.copy()
    plot_df["starttime"] = pd.to_numeric(plot_df["starttime"], errors="coerce")
    plot_df["endtime"] = pd.to_numeric(plot_df["endtime"], errors="coerce")
    plot_df = plot_df.dropna(subset=["starttime", "endtime"])
    plot_df = plot_df[plot_df["endtime"] >= plot_df["starttime"]]
    if plot_df.empty:
        return f"<h2>{escape(title)}</h2><p>Gantt chart unavailable: no valid timestamp ranges were found.</p>"

    plot_df = plot_df.sort_values(["starttime", "endtime"])
    rule_order = (
        plot_df.groupby("display_rule", sort=False)["starttime"]
        .min()
        .sort_values()
        .index.tolist()
    )
    rule_positions = {rule: idx for idx, rule in enumerate(rule_order)}
    rule_labels = _shortest_unique_rule_labels(rule_order)
    min_time = float(plot_df["starttime"].min())
    max_time = float(plot_df["endtime"].max())
    span = max(max_time - min_time, 1.0)

    width = 1240
    left_margin = 320
    right_margin = 60
    top_margin = 30
    row_height = 30
    bar_height = 18
    chart_width = width - left_margin - right_margin
    height = top_margin + len(rule_order) * row_height + 40

    svg_parts = [
        f"<h2>{escape(title)}</h2>",
        _csv_download_link(plot_df, f"{title.lower().replace(' ', '_')}.csv"),
        _plot_description_html(
            "Timeline of jobs. Each bar is one sample or group execution, positioned by its start and end time."
        ),
        f"<svg id='gantt-chart' viewBox='0 0 {width} {height}' class='plot' role='img' aria-label='{escape(title)}'>",
    ]

    # Axis track
    svg_parts.append(
        f"<line x1='{left_margin}' y1='{height - 24}' x2='{left_margin + chart_width}' y2='{height - 24}' class='axis-line'></line>"
    )
    svg_parts.append(
        f"<text x='{left_margin}' y='{height - 30}' class='plot-axis'>Run window: {format_seconds(span)}</text>"
    )
    tick_count = 6
    for tick_index in range(tick_count + 1):
        fraction = tick_index / tick_count
        tick_x = left_margin + chart_width * fraction
        tick_time = min_time + span * fraction
        tick_label = datetime.fromtimestamp(tick_time).strftime("%H:%M:%S")
        svg_parts.extend(
            [
                f"<line x1='{tick_x:.2f}' y1='{height - 24}' x2='{tick_x:.2f}' y2='{height - 18}' class='axis-line'></line>",
                f"<text x='{tick_x:.2f}' y='{height - 2}' text-anchor='middle' class='plot-axis'>{escape(tick_label)}</text>",
            ]
        )
    for rule, row_index in rule_positions.items():
        y = top_margin + row_index * row_height
        svg_parts.append(f"<g class='gantt-row' data-rule='{escape(str(rule))}'>")
        svg_parts.extend(
            [
                f"<text x='{left_margin - 10}' y='{y + 13}' text-anchor='end' class='plot-label'>{escape(rule_labels.get(str(rule), str(rule)))}</text>",
                f"<rect x='{left_margin}' y='{y}' width='{chart_width}' height='{bar_height}' rx='5' ry='5' class='plot-track'></rect>",
            ]
        )

        for row in plot_df.loc[plot_df["display_rule"] == rule].itertuples(index=False):
            label = str(row.entity)
            assay = str(row.assay)
            start = float(row.starttime)
            end = float(row.endtime)
            offset = start - min_time
            duration = max(end - start, 0.5)
            x = left_margin + chart_width * (offset / span)
            bar_width = max(chart_width * (duration / span), 2)
            label_text = escape(label)
            svg_parts.append(
                f"<rect x='{x:.2f}' y='{y}' width='{bar_width:.2f}' height='{bar_height}' rx='5' ry='5' "
                f"class='plot-bar-alt gantt-bar assay-{escape(assay)}' data-assay='{escape(assay)}'>"
                f"<title>{escape(assay)} | {escape(str(rule))} | {label_text} | {format_seconds(duration)}</title>"
                "</rect>"
            )
        svg_parts.append("</g>")

    svg_parts.append("</svg>")
    return "".join(svg_parts)


def _assay_legend_html(df: pd.DataFrame) -> str:
    """Render a compact legend for assay colors used in the Gantt chart."""
    if df.empty or "assay" not in df.columns:
        return ""

    assays = [assay for assay in df["assay"].dropna().unique().tolist() if assay]
    if not assays:
        return ""

    assays = sorted(assays)
    items = [
        (
            "<button type='button' class='legend-item legend-toggle is-active' data-assay='all'>"
            "<span class='legend-swatch legend-swatch-all'></span>"
            "<span class='legend-label'>All Assays</span>"
            "</button>"
        )
    ]
    items.extend(
        (
            "<button type='button' class='legend-item legend-toggle is-active' "
            f"data-assay='{escape(assay)}'>"
            f"{_assay_swatch_html(assay)}"
            f"<span class='legend-label'>{escape(_assay_display_name(assay))}</span>"
            "</button>"
        )
        for assay in assays
    )
    return f"<div class='legend gantt-legend'>{''.join(items)}</div>"


def compute_assay_output_sizes(output_root: Path) -> pd.DataFrame:
    """Compute total output size for each assay directory under seqnado_output."""
    if not output_root.exists() or not output_root.is_dir():
        return pd.DataFrame(columns=["assay", "output_size_bytes"])

    rows: list[dict[str, object]] = []
    for path in sorted(output_root.iterdir()):
        if not path.is_dir():
            continue
        if not (path / ".benchmark").exists():
            continue

        total_size = 0
        for child in path.rglob("*"):
            if not child.is_file():
                continue
            if ".benchmark" in child.parts:
                continue
            try:
                total_size += child.stat().st_size
            except OSError:
                continue

        rows.append({"assay": path.name, "output_size_bytes": total_size})

    if not rows:
        return pd.DataFrame(columns=["assay", "output_size_bytes"])

    return (
        pd.DataFrame(rows)
        .sort_values("output_size_bytes", ascending=False)
        .reset_index(drop=True)
    )


def _prepare_report_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    """Return a report-ready dataframe with derived metrics added."""
    prepared = df.copy()
    if "io_in" in prepared.columns or "io_out" in prepared.columns:
        prepared["io_signed"] = pd.to_numeric(
            prepared.get("io_out", 0), errors="coerce"
        ).fillna(0) - pd.to_numeric(prepared.get("io_in", 0), errors="coerce").fillna(0)
    return prepared


def _format_gb_from_mb(value_mb: float) -> str:
    """Format a megabyte-valued metric as decimal GB."""
    return f"{value_mb / 1000:,.2f} GB"


def _format_gb_from_bytes(value_bytes: float) -> str:
    """Format a byte-valued metric as decimal GB."""
    return f"{value_bytes / 1_000_000_000:,.2f} GB"


def _summary_cards(summary: BenchmarkSummary) -> list[tuple[str, str]]:
    """Return overview card label/value pairs."""
    return [
        ("Jobs", f"{summary.jobs:,d}"),
        ("Rules", f"{summary.groups:,d}"),
        ("Total Runtime", format_seconds(summary.total_runtime_seconds)),
        ("Total CPU Time", format_seconds(summary.total_cpu_time_seconds)),
        ("Peak RSS", _format_gb_from_mb(summary.peak_rss_mb)),
        ("Total I/O In", _format_gb_from_bytes(summary.total_io_in)),
        ("Total I/O Out", _format_gb_from_bytes(summary.total_io_out)),
    ]


def _summary_cards_html(summary: BenchmarkSummary) -> str:
    """Render summary cards for the overview section."""
    return "".join(
        f"<div class='card'><div class='label'>{escape(label)}</div><div class='value'>{escape(value)}</div></div>"
        for label, value in _summary_cards(summary)
    )


def _report_sections(
    df: pd.DataFrame,
    top_n: int,
    assay_sizes_plot: str,
    timeline_df: pd.DataFrame,
    read_counts_df: pd.DataFrame,
) -> list[tuple[str, str, str]]:
    """Build the ordered report sections."""
    sections = [
        (
            "run-timeline",
            "Run Timeline",
            f"{_assay_legend_html(timeline_df)}{_gantt_plot_html(timeline_df, 'Run Timeline')}",
        ),
        ("output-size", "Output Size", assay_sizes_plot),
        (
            "runtime",
            "Runtime",
            _box_plot_html(df, "s", "Runtime", "Runtime (s)", top_n=top_n),
        ),
        (
            "max-rss",
            "Max RSS",
            _box_plot_html(df, "max_rss", "Max RSS", "Max RSS (MB)", top_n=top_n),
        ),
    ]
    if _metric_available(df, "max_uss"):
        sections.append(
            (
                "max-uss",
                "Max USS",
                _box_plot_html(df, "max_uss", "Max USS", "Max USS (MB)", top_n=top_n),
            )
        )
    if _metric_available(df, "max_pss"):
        sections.append(
            (
                "max-pss",
                "Max PSS",
                _box_plot_html(df, "max_pss", "Max PSS", "Max PSS (MB)", top_n=top_n),
            )
        )
    if _metric_available(df, "cpu_usage"):
        sections.append(
            (
                "cpu-usage",
                "CPU Usage",
                _box_plot_html(
                    df, "cpu_usage", "CPU Usage", "CPU Usage (%)", top_n=top_n
                ),
            )
        )
    if _metric_available(df, "input_size_mb"):
        sections.append(
            (
                "input-size",
                "Input Size",
                _box_plot_html(
                    df, "input_size_mb", "Input Size", "Input Size (MB)", top_n=top_n
                ),
            )
        )
    if _metric_available(df, "io_signed"):
        sections.append(
            ("io", "I/O", _box_plot_html(df, "io_signed", "I/O", "I/O", top_n=top_n))
        )
    if not read_counts_df.empty:
        sections.append(
            (
                "read-counts",
                "Read Counts",
                _read_count_line_plot_html(read_counts_df),
            )
        )
    return sections


def write_html_report(
    df: pd.DataFrame,
    benchmark_dir: Path,
    output_file: Path,
    top_n: int = 20,
    assay_sizes: pd.DataFrame | None = None,
    timeline_df: pd.DataFrame | None = None,
    read_counts_df: pd.DataFrame | None = None,
) -> None:
    """Write a self-contained HTML benchmark report."""
    summary = summarize_benchmarks(df)
    df = _prepare_report_dataframe(df)
    cards_html = _summary_cards_html(summary)

    assay_sizes = (
        assay_sizes
        if assay_sizes is not None
        else pd.DataFrame(columns=["assay", "output_size_bytes"])
    )
    timeline_df = (
        timeline_df
        if timeline_df is not None
        else pd.DataFrame(columns=["label", "starttime", "endtime"])
    )
    read_counts_df = (
        read_counts_df
        if read_counts_df is not None
        else parse_alignment_processing_logs(
            benchmark_dir if benchmark_dir.name == "seqnado_output" else benchmark_dir
        )
    )
    assay_counts_html = _overview_assay_counts_html(benchmark_dir, timeline_df)
    provenance_html = _overview_provenance_html()
    assay_sizes_plot = _bar_plot_html(
        assay_sizes,
        "assay",
        "output_size_bytes",
        "Output Size",
        "Output size",
        top_n=max(len(assay_sizes), 1),
    )
    sections = _report_sections(
        df, top_n, assay_sizes_plot, timeline_df, read_counts_df
    )
    toc_html = "".join(
        f"<a class='toc-link' href='#{escape(section_id)}'>{escape(title)}</a>"
        for section_id, title, _content in sections
    )
    sections_html = "".join(
        _collapsible_section_html(section_id, title, content)
        for section_id, title, content in sections
    )

    body = "\n".join(
        [
            "<!doctype html>",
            "<html lang='en'>",
            "<head>",
            "  <meta charset='utf-8'>",
            "  <meta name='viewport' content='width=device-width, initial-scale=1'>",
            f"  <title>SeqNado Benchmark Report</title>",
            "  <style>",
            "    :root { --bg: #f7f4ee; --panel: #fffdfa; --ink: #1f2933; --muted: #6b7280; --line: #d8d1c5; --accent: #b45309; }",
            "    body { margin: 0; font-family: 'Helvetica Neue', Helvetica, Arial, sans-serif; background: linear-gradient(180deg, #f7f4ee 0%, #efe7da 100%); color: var(--ink); }",
            "    main { max-width: 1180px; margin: 0 auto; padding: 32px 20px 56px; }",
            "    h1, h2 { margin: 0 0 14px; }",
            "    h1 { font-size: clamp(1.5rem, 1.1vw + 1.25rem, 2.4rem); line-height: 1.05; letter-spacing: -0.02em; }",
            "    p.meta { color: var(--muted); margin: 6px 0 24px; }",
            "    .report-layout { display: grid; grid-template-columns: minmax(210px, 240px) minmax(0, 1fr); gap: 24px; align-items: start; }",
            "    .report-header { margin-bottom: 12px; }",
            "    .report-toc { position: sticky; top: 20px; background: rgba(255, 253, 250, 0.84); border: 1px solid var(--line); border-radius: 16px; padding: 16px; box-shadow: 0 10px 20px rgba(31, 41, 51, 0.05); backdrop-filter: blur(10px); }",
            "    .report-toc h2 { font-size: 0.95rem; margin-bottom: 12px; text-transform: uppercase; letter-spacing: 0.08em; color: var(--muted); }",
            "    .toc-link { display: block; padding: 8px 10px; margin: 0 0 6px; border-radius: 10px; color: var(--ink); text-decoration: none; font-size: 0.95rem; }",
            "    .toc-link:hover, .toc-link.is-active { background: #fbf6ec; color: var(--accent); }",
            "    .report-content { min-width: 0; }",
            "    .overview-subtitle { font-size: 0.86rem; font-weight: 700; text-transform: uppercase; letter-spacing: 0.08em; color: var(--muted); margin: 4px 0 12px; }",
            "    .overview-assay-counts { margin-top: -6px; margin-bottom: 4px; }",
            "    .overview-provenance { margin: 0 0 18px; padding: 12px 14px; background: #fbf6ec; border: 1px solid #ece6dc; border-radius: 14px; }",
            "    .provenance-line { color: var(--ink); font-size: 0.98rem; line-height: 1.45; }",
            "    .assay-count-grid { display: flex; flex-wrap: wrap; gap: 10px; }",
            "    .assay-count-chip { display: inline-flex; align-items: center; gap: 8px; background: #fbf6ec; border: 1px solid #ece6dc; border-radius: 999px; padding: 8px 12px; font-size: 0.95rem; }",
            "    .assay-count-label { color: var(--ink); }",
            "    .assay-count-value { font-weight: 700; color: var(--accent); min-width: 1.5ch; text-align: right; }",
            "    .cards { display: grid; grid-template-columns: repeat(auto-fit, minmax(160px, 1fr)); gap: 12px; margin: 0 0 28px; }",
            "    .card { background: var(--panel); border: 1px solid var(--line); border-radius: 14px; padding: 14px 16px; box-shadow: 0 10px 20px rgba(31, 41, 51, 0.05); }",
            "    .label { font-size: 12px; text-transform: uppercase; letter-spacing: 0.08em; color: var(--muted); margin-bottom: 8px; }",
            "    .value { font-size: 28px; font-weight: 700; color: var(--accent); }",
            "    .section { background: var(--panel); border: 1px solid var(--line); border-radius: 16px; padding: 20px 20px 20px 12px; margin-bottom: 18px; box-shadow: 0 10px 20px rgba(31, 41, 51, 0.05); overflow-x: auto; text-align: left; }",
            "    .section-body { margin-bottom: 18px; }",
            "    .report-section { margin-bottom: 14px; }",
            "    .report-section > summary { list-style: none; cursor: pointer; background: var(--panel); border: 1px solid var(--line); border-radius: 16px; padding: 16px 18px; box-shadow: 0 10px 20px rgba(31, 41, 51, 0.05); font-weight: 700; }",
            "    .report-section > summary::-webkit-details-marker { display: none; }",
            "    .report-section > summary::after { content: '+'; float: right; color: var(--accent); font-size: 1.1rem; line-height: 1; }",
            "    .report-section[open] > summary { border-bottom-left-radius: 12px; border-bottom-right-radius: 12px; margin-bottom: 8px; }",
            "    .report-section[open] > summary::after { content: '−'; }",
            "    .report-table { width: 100%; border-collapse: collapse; font-size: 14px; }",
            "    .report-table th, .report-table td { text-align: left; padding: 10px 12px; border-bottom: 1px solid #ece6dc; vertical-align: top; }",
            "    .report-table th { background: #fbf6ec; }",
            "    code { background: #fbf6ec; padding: 2px 6px; border-radius: 6px; }",
            "    .plot { width: 100%; height: auto; min-width: 760px; display: block; margin: 0 0 0 -8px; }",
            "    .plot-track { fill: #efe7da; }",
            "    .plot-bar { fill: #b45309; }",
            "    .plot-bar-alt { fill: #475569; }",
            "    .plot-bar-alt.assay-atac { fill: #c2410c; }",
            "    .plot-bar-alt.assay-chip { fill: #0f766e; }",
            "    .plot-bar-alt.assay-rna { fill: #2563eb; }",
            "    .plot-bar-alt.assay-meth { fill: #7c3aed; }",
            "    .plot-bar-alt.assay-snp { fill: #f6af2c; }",
            "    .plot-bar-alt.assay-cat { fill: #0891b2; }",
            "    .plot-bar-alt.assay-mcc { fill: #65a30d; }",
            "    .plot-bar-alt.assay-crispr { fill: #db2777; }",
            "    .plot-bar-alt.assay-multiomics { fill: #000000; }",
            "    .plot-label, .plot-value, .plot-axis { fill: #1f2933; font-size: 12px; font-family: 'Helvetica Neue', Helvetica, Arial, sans-serif; }",
            "    .plot-axis { fill: #6b7280; }",
            "    .axis-line { stroke: #9ca3af; stroke-width: 1.5; }",
            "    .legend { display: flex; flex-wrap: wrap; gap: 10px 14px; margin: 0 0 12px; }",
            "    .legend-item { display: inline-flex; align-items: center; gap: 8px; font-size: 13px; color: #1f2933; }",
            "    .legend-toggle { background: #fffdfa; border: 1px solid #d8d1c5; border-radius: 999px; padding: 7px 12px; cursor: pointer; transition: opacity 120ms ease, border-color 120ms ease, box-shadow 120ms ease; }",
            "    .legend-toggle:hover { border-color: #b89e74; box-shadow: 0 2px 8px rgba(31, 41, 51, 0.08); }",
            "    .legend-toggle.is-inactive { opacity: 0.45; }",
            "    .legend-toggle.is-active { opacity: 1; }",
            "    .legend-swatch { width: 14px; height: 14px; border-radius: 999px; border: 1px solid rgba(31, 41, 51, 0.15); display: inline-block; }",
            "    .legend-swatch-all { background: linear-gradient(90deg, #0891b2 0%, #2563eb 50%, #9333ea 100%); }",
            "    .chip { display: inline-flex; align-items: center; gap: 8px; padding: 6px 10px; border-radius: 999px; background: #fbf6ec; border: 1px solid #ece6dc; font-size: 13px; color: #1f2933; }",
            "    .chip-swatch { width: 12px; height: 12px; border-radius: 999px; border: 1px solid rgba(31, 41, 51, 0.18); display: inline-block; flex: 0 0 auto; }",
            "    .read-count-series, .read-count-point { transition: opacity 120ms ease, filter 120ms ease; }",
            "    .read-count-series.is-dimmed, .read-count-point.is-dimmed { opacity: 0.18; filter: grayscale(1); }",
            "    .plot-tooltip { position: fixed; z-index: 20; pointer-events: none; background: rgba(31, 41, 51, 0.94); color: #fffdfa; padding: 8px 10px; border-radius: 10px; box-shadow: 0 10px 20px rgba(31, 41, 51, 0.18); font-size: 12px; line-height: 1.35; max-width: 240px; }",
            "    .export-link { display: inline-block; margin: 0 0 12px; font-size: 13px; color: #8b5e00; text-decoration: none; border-bottom: 1px solid rgba(139, 94, 0, 0.35); }",
            "    .export-link:hover { color: #5f3f00; border-bottom-color: rgba(95, 63, 0, 0.55); }",
            "    .plot-description { margin: 0 0 12px; color: #5f6b76; font-size: 13px; line-height: 1.45; max-width: 72ch; }",
            "    .chart-legend { margin-bottom: 14px; }",
            "    .plot-row-band { fill: rgba(251, 246, 236, 0.55); }",
            "    .plot-grid-line { stroke: rgba(156, 163, 175, 0.22); stroke-width: 1; stroke-dasharray: 3 5; }",
            "    .plot-zero-line { stroke: rgba(180, 83, 9, 0.55); stroke-width: 2; }",
            "    .box-whisker { stroke: #64748b; stroke-width: 2; opacity: 0.9; }",
            "    .box-whisker-cap { stroke: #475569; stroke-width: 2.2; }",
            "    .box-rect { stroke: rgba(31, 41, 51, 0.28); stroke-width: 1.2; fill-opacity: 0.88; }",
            "    .box-median { stroke: #fffdfa; stroke-width: 2.4; }",
            "    .box-median-dot { fill: #1f2933; opacity: 0.82; }",
            "    .plot-track-line { stroke: rgba(231, 223, 209, 0.9); stroke-width: 1.5; }",
            "    .gantt-bar { transition: opacity 120ms ease, filter 120ms ease; }",
            "    .gantt-bar.is-dimmed { opacity: 0.18; filter: grayscale(1); }",
            "    @media (max-width: 960px) { .report-layout { grid-template-columns: 1fr; } .report-toc { position: static; margin-bottom: 18px; } }",
            "  </style>",
            "</head>",
            "<body>",
            "  <main>",
            "    <div class='report-header'>",
            "      <h1>SeqNado Benchmark Report</h1>",
            f"      <p class='meta'>Source: <code>{escape(str(benchmark_dir.resolve()))}</code></p>",
            "    </div>",
            "    <div class='report-layout'>",
            "      <aside class='report-toc'>",
            "        <h2>Contents</h2>",
            "        <a class='toc-link' href='#overview'>Overview</a>",
            f"        {toc_html}",
            "      </aside>",
            "      <div class='report-content'>",
            "        <section class='section' id='overview'>",
            f"          <div class='cards'>{cards_html}</div>",
            f"          {provenance_html}",
            f"          {assay_counts_html}",
            "        </section>",
            f"        {sections_html}",
            "      </div>",
            "    </div>",
            "    <script>",
            "      (() => {",
            "        const openSectionFromHash = () => {",
            "          const hash = window.location.hash ? window.location.hash.slice(1) : '';",
            "          if (!hash) return;",
            "          const target = document.getElementById(hash);",
            "          if (target && target.tagName.toLowerCase() === 'details') {",
            "            target.open = true;",
            "          }",
            "          document.querySelectorAll('.toc-link').forEach((link) => {",
            "            const targetHash = (link.getAttribute('href') || '').replace(/^#/, '');",
            "            link.classList.toggle('is-active', targetHash === hash);",
            "          });",
            "        };",
            "        const legend = document.querySelector('.gantt-legend');",
            "        const chart = document.getElementById('gantt-chart');",
            "        if (legend && chart) {",
            "          const toggles = Array.from(legend.querySelectorAll('.legend-toggle[data-assay]'));",
            "          const assayToggles = toggles.filter((node) => node.dataset.assay !== 'all');",
            "          const bars = Array.from(chart.querySelectorAll('.gantt-bar[data-assay]'));",
            "          const setState = (activeAssay) => {",
            "            const showAll = !activeAssay;",
            "            const allToggle = legend.querySelector('.legend-toggle[data-assay=\"all\"]');",
            "            if (allToggle) {",
            "              allToggle.classList.toggle('is-active', showAll);",
            "              allToggle.classList.toggle('is-inactive', !showAll);",
            "            }",
            "            assayToggles.forEach((toggle) => {",
            "              const isActive = showAll || toggle.dataset.assay === activeAssay;",
            "              toggle.classList.toggle('is-active', isActive);",
            "              toggle.classList.toggle('is-inactive', !isActive);",
            "            });",
            "            bars.forEach((bar) => {",
            "              const isDimmed = !showAll && bar.dataset.assay !== activeAssay;",
            "              bar.classList.toggle('is-dimmed', isDimmed);",
            "            });",
            "          };",
            "          toggles.forEach((toggle) => {",
            "            toggle.addEventListener('click', () => {",
            "              if (toggle.dataset.assay === 'all') {",
            "                legend.dataset.activeAssay = '';",
            "                setState('');",
            "                return;",
            "              }",
            "              const active = legend.dataset.activeAssay || '';",
            "              const next = active === toggle.dataset.assay ? '' : toggle.dataset.assay;",
            "              legend.dataset.activeAssay = next;",
            "              setState(next);",
            "            });",
            "          });",
            "          setState('');",
            "          }",
            "        const readCountLegend = document.querySelector('.read-count-legend');",
            "        const readCountChart = document.getElementById('read-count-chart');",
            "        const readCountTooltip = document.getElementById('read-count-tooltip');",
            "        if (readCountLegend && readCountChart) {",
            "          const toggles = Array.from(readCountLegend.querySelectorAll('.legend-toggle[data-sample]'));",
            "          const series = Array.from(readCountChart.querySelectorAll('.read-count-series[data-sample], .read-count-point[data-sample]'));",
            "          const points = Array.from(readCountChart.querySelectorAll('.read-count-point[data-sample]'));",
            "          const setReadCountState = (activeSample) => {",
            "            toggles.forEach((toggle) => {",
            "              const isActive = !activeSample || toggle.dataset.sample === activeSample;",
            "              toggle.classList.toggle('is-active', isActive);",
            "              toggle.classList.toggle('is-inactive', !isActive);",
            "            });",
            "            series.forEach((node) => {",
            "              const isDimmed = !!activeSample && node.dataset.sample !== activeSample;",
            "              node.classList.toggle('is-dimmed', isDimmed);",
            "            });",
            "          };",
            "          toggles.forEach((toggle) => {",
            "            toggle.addEventListener('click', () => {",
            "              const active = readCountLegend.dataset.activeSample || '';",
            "              const next = active === toggle.dataset.sample ? '' : toggle.dataset.sample;",
            "              readCountLegend.dataset.activeSample = next;",
            "              setReadCountState(next);",
            "            });",
            "          });",
            "          if (readCountTooltip) {",
            "            const showTooltip = (event) => {",
            "              const point = event.currentTarget;",
            "              readCountTooltip.hidden = false;",
            "              readCountTooltip.innerHTML = `${point.dataset.sample}<br>${point.dataset.step}<br>${point.dataset.reads} reads (${point.dataset.readsMillions}M)`;",
            "              readCountTooltip.style.left = `${event.clientX + 14}px`;",
            "              readCountTooltip.style.top = `${event.clientY + 14}px`;",
            "            };",
            "            const hideTooltip = () => {",
            "              readCountTooltip.hidden = true;",
            "            };",
            "            points.forEach((point) => {",
            "              point.addEventListener('mouseenter', showTooltip);",
            "              point.addEventListener('mousemove', showTooltip);",
            "              point.addEventListener('mouseleave', hideTooltip);",
            "            });",
            "          }",
            "          setReadCountState('');",
            "        }",
            "        document.querySelectorAll('.toc-link').forEach((link) => {",
            "          link.addEventListener('click', (event) => {",
            "            const targetId = (link.getAttribute('href') || '').replace(/^#/, '');",
            "            const target = document.getElementById(targetId);",
            "            if (!target) return;",
            "            event.preventDefault();",
            "            if (target.tagName.toLowerCase() === 'details') {",
            "              target.open = true;",
            "            }",
            "            if (window.history && window.history.replaceState) {",
            "              window.history.replaceState(null, '', `#${targetId}`);",
            "            } else {",
            "              window.location.hash = targetId;",
            "            }",
            "            target.scrollIntoView({ behavior: 'smooth', block: 'start' });",
            "            openSectionFromHash();",
            "          });",
            "        });",
            "        window.addEventListener('hashchange', openSectionFromHash);",
            "        openSectionFromHash();",
            "      })();",
            "    </script>",
            "  </main>",
            "</body>",
            "</html>",
        ]
    )

    output_file.parent.mkdir(parents=True, exist_ok=True)
    output_file.write_text(body, encoding="utf-8")
