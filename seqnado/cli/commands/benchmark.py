"""Generate reports from Snakemake benchmark files."""

from __future__ import annotations

from pathlib import Path

import typer
from loguru import logger

from seqnado.cli.app_instance import app
from seqnado.cli.benchmark_helpers import (
    compute_assay_output_sizes,
    discover_snakemake_logs,
    format_compact_number,
    load_benchmark_table,
    parse_alignment_processing_logs,
    parse_snakemake_logs_timeline,
    summarize_benchmarks,
    write_html_report,
)
from seqnado.cli.utils import _configure_logging, cli_print_table, verbose_option


def _resolve_benchmark_dir(benchmark_dir: Path) -> Path | None:
    """Resolve a benchmark directory, including common SeqNado default locations."""
    if benchmark_dir.exists():
        return benchmark_dir

    if benchmark_dir != Path(".benchmark"):
        return None

    candidates = [Path(".benchmark")]

    multiomics_root = Path("seqnado_output")
    if multiomics_root.exists():
        assay_benchmark_dirs = sorted(
            path / ".benchmark"
            for path in multiomics_root.iterdir()
            if path.is_dir() and (path / ".benchmark").exists()
        )
        if assay_benchmark_dirs:
            candidates.append(multiomics_root)
        candidates.extend(assay_benchmark_dirs)

    for candidate in candidates:
        if candidate.exists():
            return candidate

    return None


def _infer_output_root(benchmark_dir: Path) -> Path | None:
    """Infer the seqnado_output root from a resolved benchmark path."""
    if benchmark_dir.name == "seqnado_output":
        return benchmark_dir

    parts = benchmark_dir.parts
    if "seqnado_output" not in parts:
        return None

    seqnado_index = parts.index("seqnado_output")
    return Path(*parts[: seqnado_index + 1])


def _infer_run_root(benchmark_dir: Path) -> Path | None:
    """Infer the workflow run root that contains .snakemake/log."""
    output_root = _infer_output_root(benchmark_dir)
    if output_root is not None:
        return output_root.parent

    if benchmark_dir.name == ".benchmark":
        return benchmark_dir.parent

    return None


def _resolve_report_outputs(
    benchmark_dir: Path,
    output_html: Path,
    output_tsv: Path,
) -> tuple[Path, Path]:
    """Resolve default report output paths relative to seqnado_output when possible."""
    output_root = _infer_output_root(benchmark_dir)

    html_path = output_html
    if output_html == Path("benchmark_report.html") and output_root is not None:
        html_path = output_root / "seqnado_benchmark.html"

    tsv_path = output_tsv
    if output_tsv == Path("benchmark_summary.tsv") and output_root is not None:
        tsv_path = output_root / "seqnado_benchmark.tsv"

    return html_path, tsv_path


@app.command(help="Aggregate Snakemake .benchmark TSV files into a summary TSV and HTML report.")
def benchmark(
    benchmark_dir: Path = typer.Argument(
        Path(".benchmark"),
        exists=False,
        file_okay=False,
        dir_okay=True,
        readable=True,
        resolve_path=False,
        metavar="[BENCHMARK_DIR]",
        help="Directory containing benchmark TSV files or a SeqNado output root. Defaults to .benchmark.",
    ),
    output_html: Path = typer.Option(
        Path("benchmark_report.html"),
        "--html",
        "-o",
        help="Path to write the HTML report. Defaults to seqnado_output/seqnado_benchmark.html when available.",
    ),
    output_tsv: Path = typer.Option(
        Path("benchmark_summary.tsv"),
        "--tsv",
        help="Path to write the aggregated benchmark table. Defaults to seqnado_output/seqnado_benchmark.tsv when available.",
    ),
    top_n: int = typer.Option(
        20,
        "--top",
        min=1,
        help="Number of longest and highest-memory jobs to highlight in the HTML report.",
    ),
    verbose: bool = verbose_option(),
) -> None:
    """
    Generate a benchmark report from Snakemake benchmark output.

    Examples:
        seqnado benchmark
        seqnado benchmark seqnado_output
        seqnado benchmark seqnado_output/chip/.benchmark
        seqnado benchmark .benchmark --html qc/benchmark_report.html --tsv qc/benchmark_summary.tsv
    """
    _configure_logging(verbose)

    resolved_benchmark_dir = _resolve_benchmark_dir(benchmark_dir)
    if resolved_benchmark_dir is None:
        logger.error(
            f"Benchmark directory not found: {benchmark_dir}. "
            "Tried .benchmark and common SeqNado output locations."
        )
        raise typer.Exit(code=1)

    if resolved_benchmark_dir != benchmark_dir:
        logger.info(f"Using benchmark directory: {resolved_benchmark_dir}")

    table = load_benchmark_table(resolved_benchmark_dir)
    if table.empty:
        logger.error(f"No benchmark TSV files found under {resolved_benchmark_dir}")
        raise typer.Exit(code=1)

    output_html, output_tsv = _resolve_report_outputs(resolved_benchmark_dir, output_html, output_tsv)
    output_root = _infer_output_root(resolved_benchmark_dir)
    run_root = _infer_run_root(resolved_benchmark_dir)
    assay_sizes = compute_assay_output_sizes(output_root) if output_root is not None else None
    read_counts_df = parse_alignment_processing_logs(output_root) if output_root is not None else None
    log_candidates = discover_snakemake_logs(run_root) if run_root is not None else []
    timeline_df = parse_snakemake_logs_timeline(log_candidates) if log_candidates else None

    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(output_tsv, sep="\t", index=False)
    write_html_report(
        table,
        benchmark_dir=resolved_benchmark_dir,
        output_file=output_html,
        top_n=top_n,
        assay_sizes=assay_sizes,
        timeline_df=timeline_df,
        read_counts_df=read_counts_df,
    )

    summary = summarize_benchmarks(table)
    rows = [
        ["Jobs", f"{summary.jobs:,d}"],
        ["Total Runtime (min)", f"{summary.total_runtime_seconds / 60:,.2f}"],
        ["Total CPU Time (min)", f"{summary.total_cpu_time_seconds / 60:,.2f}"],
        ["Peak RSS (GB)", f"{summary.peak_rss_mb / 1000:,.2f}"],
        ["Total I/O In (GB)", f"{summary.total_io_in / 1_000_000_000:,.2f}"],
        ["Total I/O Out (GB)", f"{summary.total_io_out / 1_000_000_000:,.2f}"],
    ]
    cli_print_table(["Metric", "Value"], rows)

    logger.success(f"Wrote HTML report to {output_html}")
    logger.success(f"Wrote aggregated TSV to {output_tsv}")
