"""Helpers for aggregating Snakemake benchmark TSV files into a report."""

from __future__ import annotations

from dataclasses import dataclass
from html import escape
from pathlib import Path
import re
from datetime import datetime
from urllib.parse import quote

import pandas as pd


NUMERIC_COLUMNS = (
    "s",
    "h:m:s",
    "max_rss",
    "max_vms",
    "max_uss",
    "max_pss",
    "io_in",
    "io_out",
    "mean_load",
    "cpu_time",
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


def discover_benchmark_files(benchmark_dir: Path) -> list[Path]:
    """Return all benchmark TSV files under a benchmark directory."""
    return sorted(
        p for p in benchmark_dir.rglob("*.tsv") if p.is_file() and ".benchmark" in p.parts
    )


def discover_snakemake_logs(run_root: Path) -> list[Path]:
    """Return Snakemake log files under a workflow run root."""
    log_dir = run_root / ".snakemake" / "log"
    if not log_dir.exists():
        return []
    return sorted(p for p in log_dir.glob("*.snakemake.log") if p.is_file())


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
            columns=["jobid", "rule", "display_rule", "assay", "label", "entity", "starttime", "endtime"]
        )

    timeline = pd.DataFrame(rows)
    timeline = timeline[timeline["rule"] != "all"]
    return timeline.sort_values(["starttime", "endtime", "jobid"]).reset_index(drop=True)


def parse_snakemake_logs_timeline(log_files: list[Path]) -> pd.DataFrame:
    """Parse and merge multiple Snakemake logs, keeping the most recent rule/entity entry."""
    if not log_files:
        return pd.DataFrame(
            columns=["jobid", "rule", "display_rule", "assay", "label", "entity", "starttime", "endtime"]
        )

    timelines = [parse_snakemake_log_timeline(log_file) for log_file in log_files]
    timelines = [timeline for timeline in timelines if not timeline.empty]
    if not timelines:
        return pd.DataFrame(
            columns=["jobid", "rule", "display_rule", "assay", "label", "entity", "starttime", "endtime"]
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
        frame["group"] = str(display_path.parent) if display_path.parent != Path(".") else "root"
        frame["job"] = benchmark_file.stem
        assay, display_rule = _split_benchmark_path(display_path)
        frame["assay"] = assay
        frame["display_rule"] = display_rule
        rows.append(frame)

    if not rows:
        return pd.DataFrame()

    combined = pd.concat(rows, ignore_index=True)
    if "s" in combined.columns:
        combined = combined.sort_values(["s", "benchmark_file"], ascending=[False, True])
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
                if len(other_parts) >= suffix_len and "/".join(other_parts[-suffix_len:]) == candidate:
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
        total_runtime_seconds=float(df["s"].fillna(0).sum()) if "s" in df.columns else 0.0,
        total_cpu_time_seconds=float(df["cpu_time"].fillna(0).sum()) if "cpu_time" in df.columns else 0.0,
        peak_rss_mb=float(df["max_rss"].fillna(0).max()) if "max_rss" in df.columns else 0.0,
        total_io_in=float(df["io_in"].fillna(0).sum()) if "io_in" in df.columns else 0.0,
        total_io_out=float(df["io_out"].fillna(0).sum()) if "io_out" in df.columns else 0.0,
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
            "io_in": "IO In",
            "io_out": "IO Out",
            "mean_load": "Mean Load",
        }
    )
    for column in ("Runtime (s)", "CPU Time (s)", "Max RSS (MB)", "Mean Load"):
        if column in display.columns:
            display[column] = display[column].map(_format_metric)
    for column in ("IO In", "IO Out"):
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
        value = float(getattr(row, value_column))
        y = top_margin + idx * row_height
        bar_width = 0 if max_value == 0 else chart_width * (value / max_value)
        bar_fill = ASSAY_COLORS.get(str(label), ASSAY_COLORS["other"]) if label_column == "assay" else "#b45309"
        display_value = (
            format_compact_number(value)
            if value_column in {"io_in", "io_out"}
            else format_bytes(value)
            if value_column == "output_size_bytes"
            else _format_metric(value)
        )
        svg_parts.extend(
            [
                f"<text x='{left_margin - 10}' y='{y + 13}' text-anchor='end' class='plot-label'>{escape(label)}</text>",
                f"<rect x='{left_margin}' y='{y}' width='{chart_width}' height='{bar_height}' rx='6' ry='6' class='plot-track'></rect>",
                f"<rect x='{left_margin}' y='{y}' width='{bar_width:.2f}' height='{bar_height}' rx='6' ry='6' class='plot-bar' style='fill:{bar_fill}'><title>{escape(label)} | {escape(display_value)}</title></rect>",
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
        tick_label = (
            format_compact_number(tick_value)
            if value_column in {"io_in", "io_out"}
            else format_bytes(tick_value)
            if value_column == "output_size_bytes"
            else _format_metric(tick_value)
        )
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

    selected = plot_df.merge(group_summary[["assay", "display_rule"]], on=["assay", "display_rule"], how="inner")
    order = list(group_summary.itertuples(index=False, name=None))
    stats_rows: list[dict[str, object]] = []
    for assay, display_rule, _, _ in order:
        values = selected.loc[
            (selected["assay"] == assay) & (selected["display_rule"] == display_rule), value_column
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
        .index
        .tolist()
    )
    rule_positions = {rule: idx for idx, rule in enumerate(rule_order)}
    rule_labels = _shortest_unique_rule_labels(rule_order)
    present_assays = [assay for assay in ASSAY_COLORS if assay != "other" and assay in stats_df["assay"].unique()]
    if "other" in stats_df["assay"].unique():
        present_assays.append("other")
    assay_spacing = 10
    center_index = (len(present_assays) - 1) / 2 if present_assays else 0
    assay_offsets = {
        assay: int(round((idx - center_index) * assay_spacing))
        for idx, assay in enumerate(present_assays)
    }
    max_assays_per_rule = int(stats_df.groupby("display_rule")["assay"].nunique().max())
    max_value = float(stats_df["max"].max())
    if max_value <= 0:
        max_value = 1.0

    width = 1040
    left_margin = 300
    right_margin = 70
    top_margin = 28
    row_height = max(34, 18 + max_assays_per_rule * assay_spacing)
    box_height = 8
    bottom_margin = 46
    chart_width = width - left_margin - right_margin
    height = top_margin + len(rule_order) * row_height + bottom_margin
    tick_count = 5

    def scale_x(value: float) -> float:
        return left_margin + chart_width * (value / max_value)

    value_formatter = (
        format_compact_number
        if value_column in {"io_in", "io_out"}
        else format_bytes
        if value_column == "output_size_bytes"
        else _format_metric
    )

    svg_parts = [
        f"<h2>{escape(title)}</h2>",
        _csv_download_link(selected, f"{title.lower().replace(' ', '_')}.csv"),
        _plot_description_html(
            f"Distribution of {value_label.lower()} for the highest-signal rule and assay combinations. Each box summarizes all matching benchmark records with min, lower quartile, median, upper quartile, and max."
        ),
        f"<svg viewBox='0 0 {width} {height}' class='plot' role='img' aria-label='{escape(title)}'>",
    ]

    for rule, row_index in rule_positions.items():
        center_y = top_margin + row_index * row_height + (row_height / 2)
        svg_parts.extend(
            [
                f"<text x='{left_margin - 10}' y='{center_y + 4:.2f}' text-anchor='end' class='plot-label'>{escape(rule_labels.get(str(rule), str(rule)))}</text>",
                f"<line x1='{left_margin}' y1='{center_y:.2f}' x2='{left_margin + chart_width}' y2='{center_y:.2f}' class='plot-track-line'></line>",
            ]
        )

    for idx, row in enumerate(stats_df.itertuples(index=False), start=0):
        center_y = top_margin + rule_positions[str(row.display_rule)] * row_height + (row_height / 2)
        y = center_y - (box_height / 2) + assay_offsets.get(str(row.assay), 0)
        color = ASSAY_COLORS.get(str(row.assay), ASSAY_COLORS["other"])
        min_x = scale_x(float(row.min))
        q1_x = scale_x(float(row.q1))
        median_x = scale_x(float(row.median))
        q3_x = scale_x(float(row.q3))
        max_x = scale_x(float(row.max))
        label = f"{row.display_rule} | {ASSAY_DISPLAY_NAMES.get(str(row.assay), str(row.assay))}"
        tooltip = (
            f"{label} | n={row.count} | min={value_formatter(row.min)} | "
            f"q1={value_formatter(row.q1)} | median={value_formatter(row.median)} | "
            f"q3={value_formatter(row.q3)} | max={value_formatter(row.max)}"
        )
        svg_parts.extend(
            [
                f"<line x1='{min_x:.2f}' y1='{y + box_height / 2:.2f}' x2='{max_x:.2f}' y2='{y + box_height / 2:.2f}' class='box-whisker'></line>",
                f"<line x1='{min_x:.2f}' y1='{y + 3:.2f}' x2='{min_x:.2f}' y2='{y + box_height - 3:.2f}' class='box-whisker'></line>",
                f"<line x1='{max_x:.2f}' y1='{y + 3:.2f}' x2='{max_x:.2f}' y2='{y + box_height - 3:.2f}' class='box-whisker'></line>",
                f"<rect x='{q1_x:.2f}' y='{y:.2f}' width='{max(q3_x - q1_x, 2):.2f}' height='{box_height}' rx='4' ry='4' fill='{color}' class='box-rect'><title>{escape(tooltip)}</title></rect>",
                f"<line x1='{median_x:.2f}' y1='{y:.2f}' x2='{median_x:.2f}' y2='{y + box_height:.2f}' class='box-median'></line>",
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
        tick_label = value_formatter(tick_value)
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


def _gantt_plot_html(df: pd.DataFrame, title: str) -> str:
    """Render a Gantt-style timeline with rule rows and sample/group bars."""
    required = {"display_rule", "assay", "entity", "starttime", "endtime"}
    if df.empty or not required.issubset(df.columns):
        return f"<h2>{escape(title)}</h2><p>Gantt chart unavailable: no Snakemake run timeline could be parsed.</p>"

    plot_df = df.loc[:, ["display_rule", "assay", "entity", "starttime", "endtime"]].dropna()
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
        .index
        .tolist()
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
            "Timeline of logged jobs parsed from Snakemake run logs. Each bar is one sample or group execution, positioned by its real wall-clock start and end time."
        ),
        "<p class='meta'>Rules are shown on the y-axis, real wall-clock time is shown on the x-axis, and each sample/group execution is drawn as its own bar.</p>",
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
            f"<span class='legend-swatch' style='background:{escape(ASSAY_COLORS.get(assay, ASSAY_COLORS['other']))}'></span>"
            f"<span class='legend-label'>{escape(ASSAY_DISPLAY_NAMES.get(assay, assay))}</span>"
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

    return pd.DataFrame(rows).sort_values("output_size_bytes", ascending=False).reset_index(drop=True)


def write_html_report(
    df: pd.DataFrame,
    benchmark_dir: Path,
    output_file: Path,
    top_n: int = 20,
    assay_sizes: pd.DataFrame | None = None,
    timeline_df: pd.DataFrame | None = None,
) -> None:
    """Write a self-contained HTML benchmark report."""
    summary = summarize_benchmarks(df)

    cards = [
        ("Jobs", f"{summary.jobs:,d}"),
        ("Groups", f"{summary.groups:,d}"),
        ("Total Runtime", format_seconds(summary.total_runtime_seconds)),
        ("Total CPU Time", format_seconds(summary.total_cpu_time_seconds)),
        ("Peak RSS", f"{summary.peak_rss_mb:,.2f} MB"),
        ("Total IO In", format_compact_number(summary.total_io_in)),
        ("Total IO Out", format_compact_number(summary.total_io_out)),
    ]

    cards_html = "".join(
        f"<div class='card'><div class='label'>{escape(label)}</div><div class='value'>{escape(value)}</div></div>"
        for label, value in cards
    )

    assay_sizes = assay_sizes if assay_sizes is not None else pd.DataFrame(columns=["assay", "output_size_bytes"])
    timeline_df = timeline_df if timeline_df is not None else pd.DataFrame(columns=["label", "starttime", "endtime"])
    assay_sizes_plot = _bar_plot_html(
        assay_sizes,
        "assay",
        "output_size_bytes",
        "Output Size By Assay",
        "Output size",
        top_n=max(len(assay_sizes), 1),
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
            "    p.meta { color: var(--muted); margin: 6px 0 24px; }",
            "    .cards { display: grid; grid-template-columns: repeat(auto-fit, minmax(160px, 1fr)); gap: 12px; margin: 0 0 28px; }",
            "    .card { background: var(--panel); border: 1px solid var(--line); border-radius: 14px; padding: 14px 16px; box-shadow: 0 10px 20px rgba(31, 41, 51, 0.05); }",
            "    .label { font-size: 12px; text-transform: uppercase; letter-spacing: 0.08em; color: var(--muted); margin-bottom: 8px; }",
            "    .value { font-size: 28px; font-weight: 700; color: var(--accent); }",
            "    .section { background: var(--panel); border: 1px solid var(--line); border-radius: 16px; padding: 20px; margin-bottom: 18px; box-shadow: 0 10px 20px rgba(31, 41, 51, 0.05); overflow-x: auto; }",
            "    .report-table { width: 100%; border-collapse: collapse; font-size: 14px; }",
            "    .report-table th, .report-table td { text-align: left; padding: 10px 12px; border-bottom: 1px solid #ece6dc; vertical-align: top; }",
            "    .report-table th { background: #fbf6ec; }",
            "    code { background: #fbf6ec; padding: 2px 6px; border-radius: 6px; }",
            "    .plot { width: 100%; height: auto; min-width: 760px; }",
            "    .plot-track { fill: #efe7da; }",
            "    .plot-bar { fill: #b45309; }",
            "    .plot-bar-alt { fill: #475569; }",
            "    .plot-bar-alt.assay-atac { fill: #c2410c; }",
            "    .plot-bar-alt.assay-chip { fill: #0f766e; }",
            "    .plot-bar-alt.assay-rna { fill: #2563eb; }",
            "    .plot-bar-alt.assay-meth { fill: #7c3aed; }",
            "    .plot-bar-alt.assay-snp { fill: #b91c1c; }",
            "    .plot-bar-alt.assay-cat { fill: #0891b2; }",
            "    .plot-bar-alt.assay-mcc { fill: #65a30d; }",
            "    .plot-bar-alt.assay-crispr { fill: #db2777; }",
            "    .plot-bar-alt.assay-multiomics { fill: #9333ea; }",
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
            "    .export-link { display: inline-block; margin: 0 0 12px; font-size: 13px; color: #8b5e00; text-decoration: none; border-bottom: 1px solid rgba(139, 94, 0, 0.35); }",
            "    .export-link:hover { color: #5f3f00; border-bottom-color: rgba(95, 63, 0, 0.55); }",
            "    .plot-description { margin: 0 0 12px; color: #5f6b76; font-size: 13px; line-height: 1.45; max-width: 72ch; }",
            "    .box-whisker { stroke: #475569; stroke-width: 1.5; }",
            "    .box-median { stroke: #fffdfa; stroke-width: 2; }",
            "    .plot-track-line { stroke: #e7dfd1; stroke-width: 2; }",
            "  </style>",
            "</head>",
            "<body>",
            "  <main>",
            "    <h1>SeqNado Benchmark Report</h1>",
            f"    <p class='meta'>Source: <code>{escape(str(benchmark_dir.resolve()))}</code></p>",
            f"    <div class='cards'>{cards_html}</div>",
            f"    <section class='section'>{_assay_legend_html(timeline_df)}{_gantt_plot_html(timeline_df, 'Run Timeline')}</section>",
            f"    <section class='section'>{assay_sizes_plot}</section>",
            f"    <section class='section'>{_box_plot_html(df, 's', 'Runtime By Rule Box Plot', 'Runtime (s)', top_n=top_n)}</section>",
            f"    <section class='section'>{_box_plot_html(df, 'max_rss', 'Max RSS By Rule Box Plot', 'Max RSS (MB)', top_n=top_n)}</section>",
            f"    <section class='section'>{_box_plot_html(df, 'io_in', 'IO In By Rule Box Plot', 'IO In', top_n=top_n)}</section>",
            f"    <section class='section'>{_box_plot_html(df, 'io_out', 'IO Out By Rule Box Plot', 'IO Out', top_n=top_n)}</section>",
            "    <script>",
            "      (() => {",
            "        const legend = document.querySelector('.gantt-legend');",
            "        const chart = document.getElementById('gantt-chart');",
            "        if (!legend || !chart) return;",
            "        const toggles = Array.from(legend.querySelectorAll('.legend-toggle[data-assay]'));",
            "        const assayToggles = toggles.filter((node) => node.dataset.assay !== 'all');",
            "        const bars = Array.from(chart.querySelectorAll('.gantt-bar[data-assay]'));",
            "        const rows = Array.from(chart.querySelectorAll('.gantt-row'));",
            "        const setState = () => {",
            "          const activeAssays = new Set(",
            "            assayToggles.filter((node) => node.classList.contains('is-active')).map((node) => node.dataset.assay)",
            "          );",
            "          const showAll = activeAssays.size === 0 || activeAssays.size === assayToggles.length;",
            "          const allToggle = legend.querySelector('.legend-toggle[data-assay=\"all\"]');",
            "          if (allToggle) {",
            "            allToggle.classList.toggle('is-active', showAll);",
            "            allToggle.classList.toggle('is-inactive', !showAll);",
            "          }",
            "          bars.forEach((bar) => {",
            "            const visible = showAll || activeAssays.has(bar.dataset.assay);",
            "            bar.style.display = visible ? '' : 'none';",
            "          });",
            "          rows.forEach((row) => {",
            "            const visibleBars = row.querySelectorAll('.gantt-bar:not([style*=\"display: none\"])');",
            "            row.style.display = visibleBars.length > 0 ? '' : 'none';",
            "          });",
            "        };",
            "        toggles.forEach((toggle) => {",
            "          toggle.addEventListener('click', () => {",
            "            if (toggle.dataset.assay === 'all') {",
            "              assayToggles.forEach((node) => {",
            "                node.classList.add('is-active');",
            "                node.classList.remove('is-inactive');",
            "              });",
            "              setState();",
            "              return;",
            "            }",
            "            toggle.classList.toggle('is-active');",
            "            toggle.classList.toggle('is-inactive', !toggle.classList.contains('is-active'));",
            "            setState();",
            "          });",
            "        });",
            "        setState();",
            "      })();",
            "    </script>",
            "  </main>",
            "</body>",
            "</html>",
        ]
    )

    output_file.parent.mkdir(parents=True, exist_ok=True)
    output_file.write_text(body, encoding="utf-8")
